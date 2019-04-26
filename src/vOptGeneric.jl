# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
__precompile__()
module vOptGeneric
using Combinatorics, Suppressor

using JuMP

export vModel,
    getvOptData,
    with_optimizer,
    optimize!,
    getvalue,
    printX_E,
    getY_N,
    vSolve,
    @variable,
    @constraint,
    @addobjective
    # writeMOP,
    # parseMOP,
    
include("MOP.jl")
include("algorithms.jl")

mutable struct vOptData
    objs::Vector{Any}    #Objectives
	objSenses::Vector{JuMP.MOI.OptimizationSense}           #Objective Senses
    Y_N::Vector{Vector{Float64}}                            #Objective values for each point
    X_E::Vector{Vector{Float64}}                            #Variable values for each point
    isacopy::Bool
end

Base.copy(vd::vOptData) = vOptData(deepcopy(vd.objs), deepcopy(vd.objSenses), [], [], true)

function getvOptData(m::Model)
    !haskey(m.ext, :vOpt) && error("This model wasn't created with vOptGeneric")
    return m.ext[:vOpt]::vOptData
end

function vModel(optimizer_factory::OptimizerFactory; args...)
    m = Model(optimizer_factory ; args...)
    m.optimize_hook = vSolve
    # m.printhook = printhook
    m.ext[:vOpt] = vOptData(Vector{QuadExpr}(), #objs
                              Vector{Symbol}(), #objSenses
                              Vector{Vector{Float64}}(), #Y_N
                              Vector{Vector{Float64}}(), #X_E
                              false) #isacopy
    return m
end

function vSolve(m::Model; relax=false, method=nothing, step = 0.5, round_results = false, verbose = true, args...)

    vd = getvOptData(m)
    if vd.isacopy
        vd.objs .= copy.(vd.objs, m)
    end

    if relax != false
        @warn "linear relaxation not yet implemented"
    end
    
    if method == :epsilon
        solve_eps(m, step, round_results, verbose ; relaxation=relax, args...)
    elseif method == :dicho || method == :dichotomy
        solve_dicho(m, round_results ; relaxation=relax, args...)
    elseif method == :Chalmet || method == :chalmet
        solve_Chalmet(m, step ; relaxation=relax, args...)
    elseif method == :lex || method == :lexico
        solve_lexico(m, verbose ; relaxation=relax, args...)
    else
        @warn("use solve(m, method = :(epsilon | dichotomy | chalmet | lexico) )")
    end

    return

end

# function printhook(io::IO, m::Model)
#     vd = getvOptData(m)
#     for i = 1:length(vd.objs)
#         println(vd.objSenses[i] == :Min ? "Min " : "Max ", vd.objs[i])
#     end
#     str = JuMP.model_str(JuMP.REPLMode, m)
#     index = findfirst("Subject to", str)
#     index !== nothing && print(str[first(index):end])
# end


    macro addobjective(m, args...)
        if length(args) != 2
            error("in @addobjective: needs three arguments: model, 
                    objective sense (Max or Min) and linear expression.")
        end
        m = esc(m)
        args[1] != :Min && args[1] != :Max && error("in @addobjective: expected Max or Min for objective sense, got $(args[1])")
        sense = args[1] == :Min ? JuMP.MOI.MIN_SENSE : JuMP.MOI.MAX_SENSE
        expr = esc(args[2])
        return quote
            f = @expression($m, $expr)
            !isa(f, JuMP.GenericAffExpr) && error("in @addobjective : vOptGeneric only supports linear objectives")
            vd = $m.ext[:vOpt]
            push!(vd.objSenses, $(esc(sense)))
            # push!(vd.objs, JuMP.MOI.ScalarAffineFunction(f))
            push!(vd.objs, f)
            f
        end
    end

# Returns coefficients for the affine part of an objective
# function prepAffObjective(m, objaff::JuMP.GenericQuadExpr)
#     # Check that no coefficients are NaN/Inf
#     JuMP.assert_isfinite(objaff)
#     if !JuMP.verify_ownership(m, objaff.aff.vars)
#         error("Variable not owned by model present in objective")
#     end
#     f = zeros(m.numCols)
#     @inbounds for ind in 1:length(objaff.aff.vars)
#         f[objaff.aff.vars[ind].col] += objaff.aff.coeffs[ind]
#     end
#     return f
# end

@deprecate print_X_E printX_E
function printX_E(m::Model)
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

function getvalue(v::VariableRef, i::Int)
    vd = getvOptData(v.m)
    return vd.X_E[i][v.col]
end

function getvalue(arr::Array{VariableRef}, i::Int)

    ret = similar(arr,Float64)
    isempty(ret) && return ret

    for I in eachindex(arr)
        v = arr[I]
        value = getvalue(v, i)
        ret[I] = value
    end

    return ret
end

getvalue(x::JuMP.Containers.DenseAxisArray, i::Int) = JuMP._map(v -> getvalue(v, i), x)

function getY_N(m::Model)
    return getvOptData(m).Y_N
end


end#module


# TODO: submit PR for :

# julia> f = JuMP.MOI.ScalarAffineFunction(@expression(m, sum(x)));

# julia> set_objective_function(m, f)

# julia> set_objective_sense(m, JuMP.MOI.MAX_SENSE, f)
# ERROR: MethodError: no method matching set_objective_sense(::Model, ::MathOptInterface.OptimizationSense, ::MathOptInterface.ScalarAffineFunction{Float64})
# Closest candidates are:
#   set_objective_sense(::Model, ::MathOptInterface.OptimizationSense) at C:\Users\Gauthier\.julia\packages\JuMP\jnmGG\src\objective.jl:45
# Stacktrace:
#  [1] top-level scope at none:0