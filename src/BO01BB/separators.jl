# This file contains functions related to cutting planes.

using CPLEX, JuMP

ϵ = 0.00001
NodeLimit = 1000

"""
Given a point `x`, generate the first closure Chavatal-Gomory cut violated by `x`.

Returns :
    - true, if a valid C-G cut is found
    - inequality coefficients, `α'*x ≤ α_0`
"""
function SP_CG_separator(x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)
    # println("n=$n ; m = $m ")
    # println("x* = $x ")

    # MIP model
    M = Model(CPLEX.Optimizer) ; JuMP.set_silent( M )
    set_optimizer_attribute(M, "CPXPARAM_MIP_Limits_Nodes", NodeLimit)

    @variable(M, α[1:n+1], Int)
    @variable(M, μ[1:m] ≥ 0)

    @constraint(M, [j in 2:n+1], α[j] ≤ μ'*A[:, j-1])
    @constraint(M, α[1] >= μ'*b -1 + 0.01) # 10^-4
    # @constraint(M, α[1] <= μ'*b ) # redundant

    @constraint(M, sum(α[j] * x[j-1] for j = 2:n+1) - α[1] >= 0.01)
    @objective(M, Max, sum(α[j] * x[j-1] for j = 2:n+1) - α[1])
    optimize!(M)

    if has_values(M)    # has feasible value
        status = JuMP.termination_status( M )
        # @info "status = $status "
        obj_val = objective_value(M)
        @info "obj_val = ", obj_val
        μ_star = value.(M[:μ])
        println(" μ_star = $μ_star ")
        α_star = round.(Int, value.(M[:α]))
        println("α_star = $(value.(M[:α]))")
        return (obj_val > 0.0+ϵ, α_star, obj_val)
    end

    return (false, zeros(n+1), 0.0)
end

"""
Given a point `x`, generate a strong Chavatal-Gomory cut violated by `x`.

Returns :
    - true, if a valid C-G cut is found
    - inequality coefficients, `α'*x ≤ α_0`
"""
function SP_CG_separator2(x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)

    # MIP model
    M = Model(CPLEX.Optimizer) ; JuMP.set_silent( M )
    set_optimizer_attribute(M, "CPXPARAM_MIP_Limits_Nodes", NodeLimit)

    @variable(M, α[1:n+1], Int)
    @variable(M, μ[1:m] ≥ 0)
    @variable(M, f[1:n+1] ≥ 0)

    @constraint(M, [j in 2:n+1], f[j] == μ'*A[:, j-1] - α[j])
    @constraint(M, f[1] == μ'*b - α[1])
    
    @constraint(M, [j in 1:n+1], f[j] ≤ 1 - 0.01)
    @constraint(M, [i in 1:m], μ[i] ≤ 1 - 0.01)

    @constraint(M, sum(α[j] * x[j-1] for j = 2:n+1) - α[1] >= 0.01)
    @objective(M, Max, sum(α[j] * x[j-1] for j = 2:n+1) - α[1])
    optimize!(M)

    if has_values(M)    # has feasible value
        status = JuMP.termination_status( M )
        # @info "status = $status "
        obj_val = objective_value(M)
        @info "obj_val = ", obj_val
        μ_star = value.(M[:μ])
        println(" μ_star = $μ_star ")
        α_star = round.(Int, value.(M[:α]))
        println("α_star = $(value.(M[:α]))")
        return (obj_val > 0.0+ϵ, α_star, obj_val)
    end

    return (false, zeros(n+1), 0.0)
end


"""
Given two points `xₗ`, `xᵣ`, generate the first closure Chavatal-Gomory cut violated by both input points.

Returns :
    - true, if a valid C-G cut is found
    - inequality coefficients, `α'*x ≤ α_0`
"""
function MP_CG_separator(x_l::Vector{Float64}, x_r::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)

    # MIP model
    M = Model(CPLEX.Optimizer) ; JuMP.set_silent( M )
    set_optimizer_attribute(M, "CPXPARAM_MIP_Limits_Nodes", NodeLimit)

    @variable(M, α[1:n+1], Int)
    @variable(M, μ[1:m] ≥ 0)
    @variable(M, l)

    @constraint(M, [j in 2:n+1], α[j] ≤ μ'*A[:, j-1])
    @constraint(M, α[1] >= μ'*b -1 + 0.01)
    # @constraint(M, α[1] <= μ'*b )

    @constraint(M, l ≤ sum(α[j] * x_l[j-1] for j = 2:n+1) - α[1])
    @constraint(M, l ≤ sum(α[j] * x_r[j-1] for j = 2:n+1) - α[1])

    @constraint(M, l >= 0.01)
    @objective(M, Max, l)
    optimize!(M)

    if has_values(M)    # has feasible value
        status = JuMP.termination_status( M )
        # @info "status = $status "
        obj_val = objective_value(M)
        @info "obj_val = ", obj_val
        μ_star = value.(M[:μ])
        println(" μ_star = $μ_star ")
        α_star = round.(Int, value.(M[:α]))
        println("α_star = $(value.(M[:α]))")
        return (obj_val > 0.0+ϵ, α_star, obj_val)
    end

    return (false, zeros(n+1), 0.0)
end


"""
Given two points `xₗ`, `xᵣ`, generate a strong Chavatal-Gomory cut violated by both input points.

Returns :
    - true, if a valid C-G cut is found
    - inequality coefficients, `α'*x ≤ α_0`
"""
function MP_CG_separator2(x_l::Vector{Float64}, x_r::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)

    # MIP model
    M = Model(CPLEX.Optimizer) ; JuMP.set_silent( M )
    set_optimizer_attribute(M, "CPXPARAM_MIP_Limits_Nodes", NodeLimit)

    @variable(M, α[1:n+1], Int)
    @variable(M, μ[1:m] ≥ 0)
    @variable(M, f[1:n+1] ≥ 0)
    @variable(M, l)

    @constraint(M, [j in 2:n+1], f[j] == μ'*A[:, j-1] - α[j])
    @constraint(M, f[1] == μ'*b - α[1])
    
    @constraint(M, [j in 1:n+1], f[j] ≤ 1 - 0.01)
    @constraint(M, [i in 1:m], μ[i] ≤ 1 - 0.01)

    @constraint(M, l ≤ sum(α[j] * x_l[j-1] for j = 2:n+1) - α[1] )
    @constraint(M, l ≤ sum(α[j] * x_r[j-1] for j = 2:n+1) - α[1] )

    @constraint(M, l >= 0.01)
    @objective(M, Max, l)
    optimize!(M)

    if has_values(M)    # has feasible value
        status = JuMP.termination_status( M )
        # @info "status = $status "
        obj_val = objective_value(M)
        @info "obj_val = ", obj_val
        μ_star = value.(M[:μ])
        println(" μ_star = $μ_star ")
        α_star = round.(Int, value.(M[:α]))
        println("α_star = $(value.(M[:α]))")
        return (obj_val > 0.0+ϵ, α_star, obj_val)
    end

    return (false, zeros(n+1), 0.0)
end




function SP_KP_heurSeparator(x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)
    cuts = []

    # for each knapsack constraint in Ax≤b 
    for i=1:m 
        # each coefficient must be positive 
        has_negative_coeff = false
        for j=1:n 
            if A[i, j] < 0.0 
                has_negative_coeff = true ; break
            end
        end
        if has_negative_coeff || b[i] < 0.0 continue end

        ratio = Dict{Int64, Float64}(j => 0.0 for j=1:n)
        for j=1:n 
            if A[i, j]>0.0 ratio[j] = (1 - x[j])/A[i, j] end
        end
        sorted_ratio = sort(collect(ratio), by=x->x[2])     # increasing order

        acc = 0.0 ; C = []
        is_valid = false
        for (j, r) in sorted_ratio
            if A[i, j]>0.0
                acc += A[i, j] ; push!(C, j)
                if acc > b[i] 
                    is_valid = true ; break
                end
            end
        end

        if !is_valid continue end 

        cut = zeros( n+1 )
        cut[1] = length(C) - 1
        for j in C
            cut[j+1] = 1.0
        end
        cut = round.(Int, cut)
        push!(cuts, cut)
    end

    return cuts
end


function SP_KP_heurSeparator2(x::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)
    cuts = []

    # for each knapsack constraint in Ax≤b 
    for i=1:m 
        # each coefficient must be positive 
        has_negative_coeff = false
        for j=1:n 
            if A[i, j] < 0.0 
                has_negative_coeff = true ; break
            end
        end
        if has_negative_coeff || b[i] < 0.0 continue end


        weights = Dict{Int64, Float64}(j => 0.0 for j=1:n)
        for j=1:n 
            if x[j]>0.0 weights[j] = A[i, j] end 
        end
        sorted_weights = sort(collect(weights), by=x->x[2], rev=true)     # decreasing order

        acc = 0.0 ; C = []
        is_valid = false
        for (j, w) in sorted_weights
            if w>0.0
                acc += w; push!(C, j)
                if acc > b[i] 
                    is_valid = true ; break
                end
            end
        end

        if !is_valid continue end 

        cut = zeros( n+1 )
        cut[1] = length(C) - 1
        for j in C
            cut[j+1] = 1.0
        end
        cut = round.(Int, cut)
        push!(cuts, cut)
    end

    return cuts
end

function MP_KP_heurSeparator2(x_l::Vector{Float64}, x_r::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)
    cuts = []

    # for each knapsack constraint in Ax≤b 
    for i=1:m 
        # each coefficient must be positive 
        has_negative_coeff = false
        for j=1:n 
            if A[i, j] < 0.0 
                has_negative_coeff = true ; break
            end
        end
        if has_negative_coeff || b[i] < 0.0 continue end

        weights_l = Dict{Int64, Float64}(j => 0.0 for j=1:n)
        weights_r = Dict{Int64, Float64}(j => 0.0 for j=1:n)
        for j=1:n 
            if x_l[j]>0.0 weights_l[j] = A[i, j] end
            if x_r[j]>0.0 weights_r[j] = A[i, j] end
        end
        sorted_weights_l = sort(collect(weights_l), by=x->x[2], rev=true)     # decreasing order
        sorted_weights_r = sort(collect(weights_r), by=x->x[2], rev=true)     # decreasing order

        acc_l = 0.0 ; acc_r = 0.0 
        cut = zeros( n+1 )

        is_valid = false
        for l = 1:n
            j = sorted_weights_l[l][1] ; w = sorted_weights_l[l][2]
            if cut[j+1] == 0.0 && w>0.0
                acc_l += w ; acc_r += w
                cut[j+1] = 1.0
            end

            if sorted_weights_r[l][1] != j && cut[sorted_weights_r[l][1] + 1] == 0.0 && sorted_weights_r[l][2]>0.0
                j = sorted_weights_r[l][1]; w = sorted_weights_r[l][2]
                acc_l += w ; acc_r += w
                cut[j+1] = 1.0
            end

            if acc_l > b[i] && acc_r > b[i]
                is_valid = true ; break
            end
        end

        if !is_valid continue end 
        
        cut[1] = sum(cut) - 1
        cut = round.(Int, cut)
        push!(cuts, cut)
    end

    return cuts
end



function MP_KP_heurSeparator(x_l::Vector{Float64}, x_r::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})
    n = size(A, 2)
    m = size(A, 1)
    cuts = []

    # for each knapsack constraint in Ax≤b 
    for i=1:m 
        # each coefficient must be positive 
        has_negative_coeff = false
        for j=1:n 
            if A[i, j] < 0.0 
                has_negative_coeff = true ; break
            end
        end
        if has_negative_coeff || b[i] < 0.0 continue end

        ratio_l = Dict{Int64, Float64}(j => 0.0 for j=1:n)
        ratio_r = Dict{Int64, Float64}(j => 0.0 for j=1:n)
        for j=1:n 
            if A[i, j]>0.0 
                ratio_l[j] = (1 - x_l[j])/A[i, j] ; ratio_r[j] = (1 - x_r[j])/A[i, j]
            end
        end
        sorted_ratio_l = sort(collect(ratio_l), by=x->x[2])     # increasing order
        sorted_ratio_r = sort(collect(ratio_r), by=x->x[2])     # increasing order

        acc_l = 0.0 ; acc_r = 0.0 
        cut = zeros( n+1 )

        is_valid = false
        for l = 1:n
            j = sorted_ratio_l[l][1] ; r = sorted_ratio_l[l][2]
            if cut[j+1] == 0.0 && A[i, j]>0.0
                acc_l += A[i, j] ; acc_r += A[i, j]
                cut[j+1] = 1.0
            end

            if sorted_ratio_r[l][1] != j && cut[sorted_ratio_r[l][1] + 1] == 0.0 && A[i, sorted_ratio_r[l][1]]>0.0
                j = sorted_ratio_r[l][1] 
                acc_l += A[i, j] ; acc_r += A[i, j]
                cut[j+1] = 1.0
            end

            if acc_l > b[i] && acc_r > b[i]
                is_valid = true ; break
            end
        end

        if !is_valid continue end 
        
        cut[1] = sum(cut) - 1
        cut = round.(Int, cut)
        push!(cuts, cut)
    end

    return cuts
end