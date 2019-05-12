# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
__precompile__()
module vOptGeneric
using Combinatorics, Suppressor

import JuMP

export vModel,
    getvOptData,
    printX_E,
    getY_N,
    vSolve,
    @addobjective
    # writeMOP,
    # parseMOP,
    
# include("MOP.jl")
include("algorithms.jl")

mutable struct vOptData
    objs::Vector{JuMP.GenericAffExpr}    #Objectives
	objSenses::Vector{JuMP.MOI.OptimizationSense}           #Objective Senses
    Y_N::Vector{Vector{Float64}}                            #Objective values for each point
    X_E::Vector{Vector{Float64}}                            #Variable values for each point
end

vOptData() = vOptData([], [], [], [])

function JuMP.copy_extension_data(data::vOptData, new_model::JuMP.AbstractModel, model::JuMP.AbstractModel)
    new_objs = [JuMP.copy(obj, new_model) for obj in data.objs]
    return vOptData(new_objs, copy(data.objSenses), [], [])
end

function getvOptData(m::JuMP.Model)
    !haskey(m.ext, :vOpt) && error("This model wasn't created with vOptGeneric")
    return m.ext[:vOpt]::vOptData
end

function vModel(optimizer_factory::JuMP.OptimizerFactory; args...)
    m = JuMP.Model(optimizer_factory ; args...)
    m.optimize_hook = vSolve
    # m.printhook = printhook
    m.ext[:vOpt] = vOptData()
    return m
end

function vSolve(m::JuMP.Model, optimizer_factory::Union{Nothing, JuMP.OptimizerFactory}=nothing; relax=false, method=nothing, step = 0.5, round_results = false, verbose = true, args...)

    vd = getvOptData(m)

    if relax != false
        @warn "linear relaxation not yet implemented"
    end
    
    if method == :epsilon
        solve_eps(m, optimizer_factory, step, round_results, verbose ; relaxation=relax, args...)
    elseif method == :dicho || method == :dichotomy
        solve_dicho(m, optimizer_factory, round_results ; relaxation=relax, args...)
    elseif method == :Chalmet || method == :chalmet
        solve_Chalmet(m, optimizer_factory, step ; relaxation=relax, args...)
    elseif method == :lex || method == :lexico
        solve_lexico(m, optimizer_factory, verbose ; relaxation=relax, args...)
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
        f = JuMP.@expression($m, $expr)
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

function printX_E(m::JuMP.Model)
    vd = getvOptData(m)
    for i = 1:length(vd.Y_N)
        print(vd.Y_N[i]," : ")
        for var = JuMP.all_variables(m)
            val = JuMP.value(var, i)
            if val != 0
                if JuMP.is_binary(var) || JuMP.is_integer(var)
                    print(JuMP.name(var), "=", round(Int, val), " ")
                else
                    print(JuMP.name(var)," =", val," ")
                end
            end
        end 
        println()
   end
end

function JuMP.value(v::JuMP.VariableRef, i::Int)
    vd = getvOptData(v.model)
    return vd.X_E[i][v.index.value]
end

function JuMP.value(arr::Array{JuMP.VariableRef}, i::Int)
    ret = similar(arr,Float64)
    for j in eachindex(arr)
        ret[j] = JuMP.value(arr[j], i)
    end
    return ret
end

function value(::AbstractArray{<:JuMP.AbstractJuMPScalar}, i::Int)
    error("`JuMP.value` is not defined for collections of JuMP types. Use" *
          " Julia's broadcast syntax instead: `JuMP.value.(x, i)`.")
end

function getY_N(m::JuMP.Model)
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