__precompile__()
module vOptGeneric

importall JuMP
export MultiModel,
    getMultiData,
    solve,
    writeMOP,
    parseMOP,
    print_X_E,
    getvalue,
    getY_N,
    @addobjective,
    @variable,
    @constraint

include("MOP.jl")

type MultiData
    objs::Array{QuadExpr}               #Objectives
	objSenses::Array{Symbol}            #Objective Senses
    Y_N::Vector{Tuple{Float64,Float64}} #Objective values for each point
    X_E::Vector{Dict{Union{Variable,Array{Variable}},Union{Float64, Array{Float64}}}}#Variable values for each point, keeps the structure of user-declared variables
    X_E_raw::Vector{Vector{Float64}}    #Variable values for each point
end

function getMultiData(m::Model)
    !haskey(m.ext, :Multi) && error("This model wasn't created with vOptGeneric")
    return m.ext[:Multi]::MultiData
end

function MultiModel(;solver=JuMP.UnsetSolver())
    m = Model(solver=solver)
    m.solvehook = solvehook
    m.ext[:Multi] = MultiData(Vector{QuadExpr}(), #objs
                              Vector{Symbol}(), #objSenses
                              Vector{NTuple{N,Float64} where N}(), #Y_N
                              Vector{Any}(), #X_E
                              Vector{Vector{Float64}}()) #X_E_raw
    return m
end

function solvehook(m::Model; suppress_warnings=false, method=nothing, step = 0.5)::Symbol
    if method == nothing
        warn("use solve(m, method = :eps)")
        return :Infeasible
    end
    if method == :eps
        return solve_eps(m, step)
    end
end

function solve_eps(m::Model, ϵ::Float64)
    any(:Cont in m.colCat) && error("Epsilon method implemented for pure integer problems only")

    #Retrieve objectives and their senses from MultiData
    md = getMultiData(m)
    f1,f2 = md.objs[1],md.objs[2]
    f1Sense,f2Sense = md.objSenses[1],md.objSenses[2]

    #Set the first objective as an objective in the JuMP model
    m.obj=f1
    m.objSense = f1Sense

    #Solve with that objective
    status = solve(m, ignore_solve_hook=true)

    #If a solution exists
    if status == :Optimal

        #Declare the epsilon-constraint (RHS will be set later)
        if f2Sense == :Min        
            @constraint(m, eps, f2.aff <= 0.0)
        else
            @constraint(m, eps, f2.aff >= 0.0)
        end

        #Dict{Symbol,Any} -> variable-name <> JuMP.Variable
        if :varDict in fieldnames(JuMP.Model)
            varDict = m.varDict
        else
            varDict = filter((k,v) -> isa(v, JuMP.Variable) || isa(v, Array{JuMP.Variable}), m.objDict)
        end

        #While a solution exists
        while status == :Optimal
            #Get the score on the objectives
            f1Val = m.objVal
            f2Val = JuMP.getvalue(f2)

            #Dict{Symbol,Any} -> variable(s) <> variable(s)-value(s)
            varValueDict = Dict(v=>JuMP.getvalue(v) for (k,v) in varDict)

            #If last solution found is dominated by this one
            if length(md.Y_N) > 0
                Y_m1 = md.Y_N[end]
                R1 = f1Sense==:Min ? (<=) : (>=)
                R2 = f1Sense==:Min ? (<=) : (>=)
                if R1(f1Val,Y_m1[1]) && R2(f2Val, Y_m1[2])
                    #Remove last solution from Y_N and X_E
                    pop!(md.Y_N) ; pop!(md.X_E) ; pop!(md.X_E_raw)
                end
            end

            #Store results in MultiData
            push!(md.Y_N, (f1Val, f2Val))
            push!(md.X_E, varValueDict)
            push!(md.X_E_raw, copy(m.colVal))

            print("z1 = ", f1Val, ", z2 = ", f2Val)
            #Set the RHS of the epsilon-constraint
            if f2Sense == :Min
                JuMP.setRHS(eps, f2Val - f2.aff.constant - ϵ)
                println(". Solving with f2 <= ", f2Val - f2.aff.constant - ϵ)
            else
                JuMP.setRHS(eps, f2Val - f2.aff.constant + ϵ)
                println(". Solving with f2 >= ", f2Val - f2.aff.constant + ϵ)
            end

            # for (k,v) in varDict
            #     setvalue(v, getvalue(v))
            # end

            #And solve again
            status = solve(m, ignore_solve_hook=true, suppress_warnings=true)
        end
        ###
        #JuMP doesn't support removing constraints
        #To leave the model unaltered we change the value of the RHS
        ###
        if f2Sense == :Min
            JuMP.setRHS(eps, Float64(typemax(Int)))
        else
            JuMP.setRHS(eps, Float64(typemin(Int)))
        end
    else
        return status
    end
    return :Optimal
end

macro addobjective(m, args...)
    if length(args) != 2
        error("in @addobjective: needs three arguments: model, 
                objective sense (Max or Min) and linear expression.")
    end
    m = esc(m)
    sense = args[1]
    if sense == :Min || sense == :Max
        sense = Expr(:quote,sense)
    end
    expr = esc(args[2])
    code = quote
        f = @expression($m, $expr)
        !isa(f, JuMP.GenericAffExpr) && error("in @addobjective : vOptGeneric only supports linear objectives")
        md = $m.ext[:Multi]
        push!(md.objSenses, $(esc(sense)))
        push!(md.objs, QuadExpr(f))
    end
    return code
end



# Returns coefficients for the affine part of an objective
function prepAffObjective(m, objaff::JuMP.GenericQuadExpr)
    # Check that no coefficients are NaN/Inf
    JuMP.assert_isfinite(objaff)
    if !JuMP.verify_ownership(m, objaff.aff.vars)
        error("Variable not owned by model present in objective")
    end
    f = zeros(m.numCols)
    @inbounds for ind in 1:length(objaff.aff.vars)
        f[objaff.aff.vars[ind].col] += objaff.aff.coeffs[ind]
    end
    return f
end

function print_X_E(m::Model)
    md = getMultiData(m)
    for i = 1:length(md.Y_N)
        print(md.Y_N[i]," : ")
        for j = 1:m.numCols
            if md.X_E_raw[i][j] != 0
                if m.colCat[j] == :Int || m.colCat[j] == :Bin
                    print(getname(m,j),"=",round(Int,md.X_E_raw[i][j])," ")
                else
                    print(getname(m,j),"=",md.X_E_raw[i][j]," ")
                end
            end
        end 
        println()
   end
end

function getvalue(arr::Array{Variable}, i::Int)
    m = first(arr).m
    md = getMultiData(m)
    return md.X_E[i][arr]
end

function getvalue(v::Variable, i::Int)
    m = v.m
    md = getMultiData(m)
    return md.X_E[i][v]
end

function getY_N(m::Model)
    return getMultiData(m).Y_N
end

end#module