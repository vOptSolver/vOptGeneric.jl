__precompile__()
module vOptGeneric

importall JuMP
import MathProgBase

export vModel,
    getvOptData,
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

type vOptData
    objs::Array{QuadExpr}               #Objectives
	objSenses::Array{Symbol}            #Objective Senses
    Y_N::Vector{NTuple{N,Float64} where N}#Objective values for each point
    X_E::Vector{Vector{Float64}}        #Variable values for each point
end

function getvOptData(m::Model)
    !haskey(m.ext, :vOpt) && error("This model wasn't created with vOptGeneric")
    return m.ext[:vOpt]::vOptData
end

function vModel(;solver=JuMP.UnsetSolver())
    m = Model(solver=solver)
    m.solvehook = solvehook
    m.ext[:vOpt] = vOptData(Vector{QuadExpr}(), #objs
                              Vector{Symbol}(), #objSenses
                              Vector{NTuple{N,Float64} where N}(), #Y_N
                              Vector{Vector{Float64}}()) #X_E
    return m
end

function solvehook(m::Model; suppress_warnings=false, method=nothing, step = 0.5)::Symbol
    if method == nothing
        warn("use solve(m, method = :(epsilon | dichotomy) )")
        return :Infeasible
    elseif method == :epsilon
        return solve_eps(m, step)
    elseif method == :dicho || method == :dichotomy
        return solve_dicho(m)
    end
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
    return quote
        f = @expression($m, $expr)
        !isa(f, JuMP.GenericAffExpr) && error("in @addobjective : vOptGeneric only supports linear objectives")
        vd = $m.ext[:vOpt]
        push!(vd.objSenses, $(esc(sense)))
        push!(vd.objs, QuadExpr(f))
    end
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
    vd = getvOptData(m)
    for i = 1:length(vd.Y_N)
        print(vd.Y_N[i]," : ")
        for j = 1:m.numCols
            if vd.X_E[i][j] != 0
                if m.colCat[j] == :Int || m.colCat[j] == :Bin
                    print(getname(m,j),"=",round(Int,vd.X_E[i][j])," ")
                else
                    print(getname(m,j),"=",vd.X_E[i][j]," ")
                end
            end
        end 
        println()
   end
end

function getvalue(v::Variable, i::Int)
    vd = getvOptData(v.m)
    return vd.X_E[i][v.col]
end

function getvalue(arr::Array{Variable}, i::Int)

    ret = similar(arr,Float64)
    isempty(ret) && return ret

    for I in eachindex(arr)
        v = arr[I]
        value = getvalue(v, i)
        ret[I] = value
    end

    return ret
end

getvalue(x::JuMP.JuMPContainer, i::Int) = JuMP._map(v -> getvalue(v, i), x)

function getY_N(m::Model)
    return getvOptData(m).Y_N
end

function solve_eps(m::Model, ϵ::Float64)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    #Set the first objective as an objective in the JuMP model
    m.obj=f1
    m.objSense = f1Sense

    #Solve with that objective
    status = solve(m, ignore_solve_hook=true)

    #If a solution exists
    if status == :Optimal

        #Declare the epsilon-constraint (RHS will be set later)
        if f2Sense == :Min        
            eps = @constraint(m, f2.aff <= 0.0)
        else
            eps = @constraint(m, f2.aff >= 0.0)
        end

        varArray = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

        #While a solution exists
        while status == :Optimal
            #Get the score on the objectives
            f1Val = m.objVal
            f2Val = JuMP.getvalue(f2)

            #If last solution found is dominated by this one
            if length(vd.Y_N) > 0
                Y_m1 = vd.Y_N[end]
                R1 = f1Sense==:Min ? (<=) : (>=)
                R2 = f1Sense==:Min ? (<=) : (>=)
                if R1(f1Val,Y_m1[1]) && R2(f2Val, Y_m1[2])
                    #Remove last solution from Y_N and X_E
                    pop!(vd.Y_N) ; pop!(vd.X_E)
                end
            end

            #Store results in vOptData
            push!(vd.Y_N, (f1Val, f2Val))
            push!(vd.X_E, JuMP.getvalue.(varArray))

            print("z1 = ", f1Val, ", z2 = ", f2Val)
            #Set the RHS of the epsilon-constraint
            if f2Sense == :Min
                JuMP.setRHS(eps, f2Val - f2.aff.constant - ϵ)
                println(". Solving with f2 <= ", f2Val - ϵ)
            else
                JuMP.setRHS(eps, f2Val - f2.aff.constant + ϵ)
                println(". Solving with f2 >= ", f2Val + ϵ)
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

function solve_dicho(m::Model)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    #Set the first objective as an objective in the JuMP model
    m.obj=f1
    m.objSense = f1Sense

    #Solve with that objective
    status = solve(m, ignore_solve_hook=true)

    #If a solution exists
    if status == :Optimal

        yr_1 = m.objVal
        yr_2 = JuMP.getvalue(f2)

        #Store results in vOptData
        push!(vd.Y_N, (yr_1, yr_2))
        push!(vd.X_E, JuMP.getvalue.(varArray))


        #Set the second objective as an objective in the JuMP model
        m.obj=f2
        m.objSense = f2Sense

        #Solve with that objective
        status = solve(m, ignore_solve_hook=true)

        if status == :Optimal

            ys_1 = JuMP.getvalue(f1)
            ys_2 = m.objVal

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, (ys_1, ys_2))
                push!(vd.X_E, JuMP.getvalue.(varArray))
                solveRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray)
            end
        
            #Lazy sort X_E and Y_N
            Y_N, X_E = vd.Y_N, vd.X_E
            s = sortperm(Y_N, by = e -> e[1])
            Y_N = Y_N[s]
            X_E = X_E[s]

            #We can use weak dominance definition here
            dom_min_min(a,b) = a[1] <= b[1] && a[2] <= b[2]
            dom_max_max(a,b) = a[1] >= b[1] && a[2] >= b[2]
            dom_min_max(a,b) = a[1] <= b[1] && a[2] >= b[2]
            dom_max_min(a,b) = a[1] >= b[1] && a[2] <= b[2]

            if f1Sense==f2Sense==:Min
                dominates = dom_min_min
            elseif f1Sense==f2Sense==:Max
                dominates = dom_max_max
            elseif f1Sense==:Min
                dominates = dom_min_max
            else
                dominates = dom_max_min
            end

            #Filter X_E and Y_N :
            inds = Int[]
            for i = 1:length(Y_N)-1
                if dominates(Y_N[i], Y_N[i+1])
                    push!(inds, i+1)
                elseif dominates(Y_N[i+1], Y_N[i])
                    push!(inds, i)
                end
            end
            deleteat!(Y_N, inds)
            deleteat!(X_E, inds)
        end
    end

    status
end

function solveRecursion(m::Model, yr_1, yr_2, ys_1, ys_2, varArray)

    vd = getvOptData(m)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    m.objSense = f1Sense

    if f1Sense==f2Sense
        lb = λ1*yr_1 + λ2*yr_2
        m.obj = λ1*f1 + λ2*f2
    else
        lb = λ1*yr_1 - λ2*yr_2
        m.obj = λ1*f1 - λ2*f2
    end

    solve(m, ignore_solve_hook=true)

    yt_1 = JuMP.getvalue(f1)
    yt_2 = JuMP.getvalue(f2)


    if f1Sense==f2Sense
        val = λ1*yt_1 + λ2*yt_2
    else
        val = λ1*yt_1 - λ2*yt_2
    end

    if f1Sense == :Min && val < lb - 1e-4 || val > lb + 1e-4
        push!(vd.Y_N, (yt_1, yt_2))
        push!(vd.X_E, JuMP.getvalue.(varArray))
        solveRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray)
        solveRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray)
    end

end

end#module