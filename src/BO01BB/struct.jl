"""
    This file contains objects types used for BO01BB algorithm.
    Consider the following multi-objective program :

    min     f1(⋅), f2(⋅)
    s.t.      Ax ≤ b
              x∈{0,1}
"""

@enum PrunedType INFEASIBILITY OPTIMALITY DOMINANCE LEAVE NONE

TOL = 10^(-6)

"""
Storage the statistics information of the BO01BB algorithm.
"""
mutable struct StatInfo
    total_times::Float64
    nb_nodes::Int64
    nb_nodes_pruned::Int64
    Gap::Float64
end

function StatInfo()
    return StatInfo(0, 0, 0, 0)
end


"""
Some parameters used in the B&B algorithm.
"""
mutable struct BBparam
    time_limit::Float64     # time limit for B&B algo
    tol::Float64            # numerical precision
    traverse::Symbol        # traverse strategy such as dfs, bfs...
    branching::Symbol       # branching strategy
end

function BBparam()
    return BBparam(300, TOL, :bfs, :arbitrary)
end


"""
Storage the components definning a bi-objective 0-1 linear program.
"""
mutable struct BO01Problem
    vd::vOptData
    varArray::Array{VariableRef}
    m::JuMP.Model
end


"""
A Sol object consists of solution x in decision space and point y in objective space.
"""
mutable struct Sol
    x::Vector{Float64}
    y::Vector{Float64}
    is_integer::Bool
end

function Sol()
    return Sol(Vector{Float64}(), Vector{Float64}(), false)
end

"""
Overload operators for the dominance order between two solutions.
"""
function Base.:<=(a::Sol, b::Sol)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1]<=b.y[1] && a.y[2]<=b.y[2]
end

function Base.:>=(a::Sol, b::Sol)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1]>=b.y[1] && a.y[2]>=b.y[2]
end

function Base.:(==)(a::Sol, b::Sol)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1]==b.y[1] && a.y[2]==b.y[2]
end

function Base.:(!=)(a::Sol, b::Sol)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1]!=b.y[1] && a.y[2]!=b.y[2]
end


"""
A vector of solutions in the natrual order (from left to right in bi-objective space).
"""
mutable struct NatrualOrderVector
    sol_vect::Vector{Sol}
end

function NatrualOrderVector()
    return NatrualOrderVector(Vector{Sol}())
end

"""
Push a solution into a vector of natrual ordered solutions, return true if it is successfully added;
or false, if it is weakly dominated by one (or more) solution of the vector. 

In case of successfully added, delete the old solutions that are weakly dominated by the new one, if have any.
"""
function Base.push!(sols::NatrualOrderVector, s::Sol)
    # add s directly if sols is empty
    if length(sols.sol_vect) == 0
        push!(sols.sol_vect, s)
        return true
    end

    # find the right position (ind) to insert, keeping the natrual order
    ind = 1
    tail = length(sols.sol_vect)
    for i in 1:tail
        # compare the first objective
        if s.y[1] < sols.sol_vect[i].y[1]
            ind += 1; continue
        else if s.y[1] == sols.sol_vect[i].y[1]
           # compare the second objective
           if s.y[2] < sols.sol_vect[i].y[2]
                break # s weakly dominates sol_vect[i]
            else
                return false # s is weakly dominated by sol_vect[i]
           end
        else
            break
        end
    end

    # insert s at position "ind"
    if tail == ind
        push!(sols.sol_vect, s)
        return true
    end

    # if s weakly dominates sol_vect[ind], replace it
    if s <= sols.sol_vect[ind]
        sols.solutions[ind] = s
    else
    # otherwise, insert s at position "ind"
        push!(sols.sol_vect, s)
        sols.solutions[ind+1:tail+1] = sols.solutions[ind:tail]
        sols.solutions[ind] = s
        return true
    end

    # find solutions (e.g. sol_vect[ind+1:supp]) weakly dominated by s (i.e. sol_vect[ind])
    supp = ind
    while supp < tail
        # if sol_vect[ind] weakly dominates sol_vect[supp+1]
        if sols.sol_vect[ind] <= sols.sol_vect[supp+1]
            supp += 1
        else
            break
        end
    end

    if supp == tail
        sols.sol_vect = sols.sol_vect[1:ind]
    else if supp > ind
        tmp = sols.sol_vect[supp+1:tail]
        sols.sol_vect = vcat(sols.sol_vect[1:ind], sols.sol_vect[supp+1:tail])
    end

    return true
end


"""
The relaxed bound set consists of segements and points.
"""
mutable struct RelaxedBoundSet
    sols::NatrualOrderVector
    segements::Vector{Vector{Int64}} # TODO : find a better structure
end


"""
The incumbent set consists in feasible solutions/points.
"""
mutable struct IncumbentSet
    sols::NatrualOrderVector
end


"""
Definition of the node object in B&B tree.
"""
mutable struct Node
    id::Int64                   # identify number (count)
    pred::Int64                 # predecessor's indices
    succs::Vector{Int64}        # successors' indices
    var::Int64                  # indice of the chosen variable to be split
    var_bound::Float64          # variable bound
    RBS::RelaxedBoundSet        # local relaxed bound set
    pareto_sols::NatrualOrderVector    # pareto optimal solutions
    actived::Bool               # if the node is active
    isPruned::Bool              # if the node is pruned
    prunedType::PrunedType      # if the node is fathomed, restore pruned type
end

# function Node()
#     # return Node(0, 0, Vector{Int64}(), 0, 0., RelaxedBoundSet(), NatrualOrderVector(), true, false, NONE)
# end

function Base.empty!(node::Node)
    
end

function prune(node::Node)
    
end

