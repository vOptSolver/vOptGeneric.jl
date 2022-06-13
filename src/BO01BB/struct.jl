"""
    This file contains objects types used for BO01BB algorithm.
    Consider the following multi-objective program :

    min     f1(⋅), f2(⋅)
    s.t.      Ax ≤ b
              x∈{0,1}
"""

@enum PrunedType NONE INFEASIBILITY OPTIMALITY DOMINANCE

TOL = 10^(-4)

"""
Storage the statistics information of the BO01BB algorithm.
"""
mutable struct StatInfo
    total_times::Float64
    nb_nodes::Int64
    nb_nodes_pruned::Int64
    Gap::Float64
    relaxation_time::Float64
    test_dom_time::Float64
    update_incumb_time::Float64
    # status::MOI.TerminationStatusCode 
end

function StatInfo()
    return StatInfo(0.0, 0, 0, 0.0, 0.0, 0.0, 0.0)
end

function Base.:show(io::IO, info::StatInfo)
    println(io, " # informations of B&B algorithm : \n",
        "total_times_used = $(info.total_times) \n",
        "total_nodes = $(info.nb_nodes) \n",
        "pruned_nodes = $(info.nb_nodes_pruned) \n",
        "GAP = $(info.Gap) \n",
        "relaxation_time = $(info.relaxation_time) \n",
        "test_dominance_time = $(info.test_dom_time) \n",
        "update_incumbent_time = $(info.update_incumb_time) ")
end


"""
Some parameters used in the B&B algorithm.
"""
mutable struct BBparam
    time_limit::Float64     # time limit for B&B algo
    traverse::Symbol        # traverse strategy such as dfs, bfs...
    branching::Symbol       # branching strategy
end

function BBparam()
    return BBparam(300, :bfs, :arbitrary)
end


"""
Storage the components definning a bi-objective 0-1 linear program.
"""
mutable struct BO01Problem
    varArray::Array{JuMP.VariableRef}
    m::JuMP.Model
    param::BBparam
    info::StatInfo
end


"""
A point `y` in the criteria space may be projected from a set of equivalent solutions `x`.
"""
mutable struct Solution
    xEquiv::Vector{Vector{Float64}}            # a set of equivalent solutions x defining the same y
    y::Vector{Float64}
    is_binary::Bool
end

function Solution()
    return Solution(Vector{Vector{Float64}}(), Vector{Float64}(), false)
end

function Solution(x::Vector{Float64}, y::Vector{Float64})
    is_binary = true
    for i = 1:length(x)
        if !(abs(x[i]-0.0) ≤ TOL || abs(x[i]-1.0) ≤ TOL)
            is_binary = false; break
        end
    end
    return Solution([x], y, is_binary)
end

"""
Add an equivalent solution associated to point `y`. 
"""
function addEquivX(sol::Solution, x::Vector{Float64})
    @assert length(x) > 0 "x cannot be empty"

    push!(sol.xEquiv, x)
    # check if x is approximately binary
    if !sol.is_binary
        for i in 1:length(x)
            if !(abs(x[i]-0.0) ≤ TOL || abs(x[i]-1.0) ≤ TOL)
                return
            end
        end
        sol.is_binary = true
    end
end

function addEquivX(sol::Solution, vecX::Vector{Vector{Float64}})
    for x in vecX
        for added_x in sol.xEquiv
            if x == added_x return end
        end
        addEquivX(sol, x)
    end
end

"""
Overload operators for the dominance order between two solutions.
"""
function Base.:show(io::IO, s::Solution)
    println(io, "Solution( \n |xEquiv| = ", length(s.xEquiv),
    "\n y = ", s.y,
    "\n is_binary ? ", s.is_binary, " )")
end


function Base.:<=(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1] ≤ b.y[1] && a.y[2] ≤ b.y[2]
end

function Base.:<(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1] < b.y[1] && a.y[2] < b.y[2]
end

function Base.:>=(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1] ≥ b.y[1] && a.y[2] ≥ b.y[2]
end

function Base.:>(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1]>b.y[1] && a.y[2]>b.y[2]
end

function Base.:(==)(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return abs(a.y[1] -b.y[1]) ≤ TOL && abs(a.y[2] - b.y[2]) ≤ TOL
end

function Base.:(!=)(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return abs(a.y[1] - b.y[1]) > TOL || abs(a.y[2] - b.y[2]) > TOL
end

function Base.isequal(a::Solution, b::Solution) # for hash table
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a == b
end


"""
Return `true` if solution `a` dominates sobution `b`; `false` otherwise.
"""
function dominate(a::Solution, b::Solution)
    return a ≤ b && a ≠ b
end


"""
A vector of solutions in the natrual order (from left to right in bi-objective space). 
"""
mutable struct NaturalOrderVector
    sols::Vector{Solution}
end

function NaturalOrderVector()
    return NaturalOrderVector(Vector{Solution}())
end

function Base.length(v::NaturalOrderVector)
    return length(v.sols)
end

function Base.:show(io::IO, nov::NaturalOrderVector)
    println(io, "NaturalOrderVector[")
    for sol in nov.sols
        print(io, sol)
    end
    println("]")
end


"""
Push a solution into a vector of natrual ordered solutions, return `true` if it is successfully added;
or `false`, if it is weakly dominated by one (or more) solution(s) in the vector. 

In case of successfully added and `filtered=true` (by defaut false), delete the old solutions that are weakly dominated by the new one.
"""
function Base.push!(natural_sols::NaturalOrderVector, sol::Solution; filtered::Bool=false)
    # add s directly if sols is empty
    if length(natural_sols) == 0
        push!(natural_sols.sols, sol)
        return true
    end

    # a binary/dichotomy search finds the location to insert 
    l = 1; r = length(natural_sols); m = 0
    while l ≤ r
        m = Int(floor((l+r)/2))
        # compare the first objective
        if sol.y[1] < natural_sols.sols[m].y[1]
            l = m+1
        elseif sol.y[1] > natural_sols.sols[m].y[1]
            r = m-1
        # in case of the equality on the first objective, compare the second obj
        elseif sol.y[2] > natural_sols.sols[m].y[2]
            l = m+1
        elseif sol.y[2] < natural_sols.sols[m].y[2]
            r  = m-1
        # in case of equality
        else
            addEquivX(natural_sols.sols[m], sol.xEquiv)
            return true
        end
    end

    if r==0 # insert at the top
        natural_sols.sols = vcat([sol], natural_sols.sols)
        m = 1
    elseif l==length(natural_sols)+1 # insert at the bottom
        push!(natural_sols.sols, sol)
        m = l
    else # insert at m
        m = m > Int(floor((l+r)/2)) ? m : m+1
        natural_sols.sols = vcat(vcat(natural_sols.sols[1:m-1], sol), natural_sols.sols[m:end])
    end

    # find points weakly dominated by the new point and delete it/them
    if filtered
        inds = []
        for i = 1:length(natural_sols)
            if dominate(sol, natural_sols.sols[i])
                push!(inds, i)
            elseif dominate(natural_sols.sols[i], sol)
                deleteat!(natural_sols.sols, m)
                return false
            end
        end
        deleteat!(natural_sols.sols, inds)
    end
    
    return true
end


"""
The relaxed bound set consists of segments and natural ordered solutions. 

Since all solutions are natural ordered, segments is defined by a dictionary where each point solution s is associated
with a boolean if s and the next/adjacent point in right forms a segment.
"""
mutable struct RelaxedBoundSet
    natural_order_vect::NaturalOrderVector
    segments::Dict{Solution, Bool}        
end

function RelaxedBoundSet()
    return RelaxedBoundSet(NaturalOrderVector(), Dict{Solution, Bool}())
end


"""
The incumbent set consists in feasible solutions.
"""
mutable struct IncumbentSet
    natural_order_vect::NaturalOrderVector
end

function IncumbentSet()
    return IncumbentSet(NaturalOrderVector())
end