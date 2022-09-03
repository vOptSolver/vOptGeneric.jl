# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
__precompile__()
module vOptGeneric
import JuMP, MathOptInterface, Combinatorics
const MOI = MathOptInterface

export vModel,
    getvOptData,
    printX_E,
    getY_N,
    vSolve,
    @addobjective,
    parseMOP,
    writeMOP

include("MOP.jl")
include("algorithms.jl")

mutable struct vOptData
    objs::Vector{JuMP.GenericAffExpr}          #Objectives
	objSenses::Vector{MOI.OptimizationSense}   #Objective Senses
    Y_N::Vector{Vector{Float64}}               #Objective values for each point
    X_E::Vector{Vector{Float64}}               #Variable values for each point
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

function vModel(optimizer_factory; args...)
    m = JuMP.Model(optimizer_factory ; args...)
    m.optimize_hook = vSolve
    # m.printhook = printhook
    m.ext[:vOpt] = vOptData()
    return m
end

function vModel(; args...)
    m = JuMP.Model(; args...)
    m.optimize_hook = vSolve
    # m.printhook = printhook
    m.ext[:vOpt] = vOptData()
    return m
end

function vSolve(m::JuMP.Model ; relax=false, method=nothing, verbose = true, kwargs...)

    vd = getvOptData(m)

    if relax != false
        @warn "linear relaxation not yet implemented"
    end

    methods = Dict{Symbol, Function}([
        :epsilon => solve_eps,
        :dicho => solve_dicho,
        :dichotomy => solve_dicho,
        :chalmet => solve_Chalmet,
        :Chalmet => solve_Chalmet,
        :lexico => solve_lexico,
        :lex => solve_lexico,
    ])

    if method in keys(methods)
        solve_function = methods[method]
        solve_function(m; verbose=verbose, relaxation=relax, kwargs...)
    else
        @warn("use solve(m, method = :(epsilon | dichotomy | chalmet | lexico))")
    end

    return m
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
    sense = args[1] == :Min ? MOI.MIN_SENSE : MOI.MAX_SENSE
    expr = esc(args[2])
    return quote
        f = JuMP.AffExpr() + JuMP.@expression($m, $expr)
        !isa(f, JuMP.GenericAffExpr) && error("in @addobjective : vOptGeneric only supports linear objectives")
        vd = $m.ext[:vOpt]
        push!(vd.objSenses, $(esc(sense)))
        push!(vd.objs, f)
        f
    end
end

function printX_E(m::JuMP.Model)
    vd = getvOptData(m)
    for i = 1:length(vd.Y_N)
        print(vd.Y_N[i]," : ")
        for var in JuMP.all_variables(m)
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

function getY_N(m::JuMP.Model)
    return getvOptData(m).Y_N
end

end
