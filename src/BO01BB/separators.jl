# This file contains functions related to cutting planes.

using CPLEX, JuMP

# include("BBtree.jl")
# include("../algorithms.jl")

ϵ = 0.00001
NodeLimit = 1000

"""
Given a point `x`, generate the strong Chavatal-Gomory cut violated by `x`.

Returns :
    - true, if a valid C-G cut is found
    - inequality coefficients, `α'*x ≤ α_0`
"""
function SP_CG_separator(x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)

    # MIP model
    M = Model(CPLEX.Optimizer) 
    # JuMP.set_silent( M )
    set_optimizer_attribute(M, "CPXPARAM_MIP_Limits_Nodes", NodeLimit)

    @variable(M, α[1:n+1], Int)
    @variable(M, μ[1:m] ≥ 0)

    @constraint(M, [j in 2:n+1], α[j] ≤ μ'*A[:, j])
    @constraint(M, α[1] > μ'*b - 1)

    @objective(M, Max, sum(α[j] * x[j] for j = 2:n+1) - α[1])
    optimize!(M)

    if has_values(M)    # has feasible value
        obj_val = objective_value(M)
        μ_star = value.(M[:μ])
        α_star = round.(Int, value.(M[:α]))
        return (obj_val > 0.0-ϵ, α_star)
    end

    return (false, zeros(n+1))
end

function SP_CG_separator2(x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)

    # MIP model
    M = Model(CPLEX.Optimizer) 
    # JuMP.set_silent( M )
    set_optimizer_attribute(M, "CPXPARAM_MIP_Limits_Nodes", NodeLimit)

    @variable(M, α[1:n+1], Int)
    @variable(M, μ[1:m] ≥ 0)
    @variable(M, f[1:n+1] ≥ 0)

    # strong valid cut
    for j = 2:n+1
        if x[j] ≈ 0.0
            @constraint(M, α[j] == 0)
        end
    end

    @constraint(M, [j in 2:n+1], f[j] == μ'*A[:, j] - α[j])
    @constraint(M, f[1] == μ'*b - α[1])
    
    @constraint(M, [j in 1:n+1], f[j] ≤ 1 - 0.01)
    @constraint(M, [i in 1:m], μ[i] ≤ 1 - 0.01)

    @objective(M, Max, sum(α[j] * x[j] for j = 2:n+1) - α[1])
    optimize!(M)

    if has_values(M)    # has feasible value
        obj_val = objective_value(M)
        μ_star = value.(M[:μ])
        α_star = round.(Int, value.(M[:α]))
        return (obj_val > 0.0-ϵ, α_star)
    end

    return (false, zeros(n+1))
end