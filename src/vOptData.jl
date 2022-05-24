

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
