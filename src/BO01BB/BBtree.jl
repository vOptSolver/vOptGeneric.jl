using JuMP

include("struct.jl")
include("cutPool.jl")

"""
Definition of the node object in B&B tree.
"""
mutable struct Node
    num::Int64                    
    depth::Int64                # depth in tree
    pred::Node                  # predecessor
    succs::Vector{Node}         # successors
    var::Int64                  # index of the chosen variable to be split
    var_bound::Int64            # variable bound
    RBS::RelaxedBoundSet        # local relaxed bound set               # TODO : erase memory
    activated::Bool             # if the node is active
    pruned::Bool                # if the node is pruned
    prunedType::PrunedType      # if the node is fathomed, restore pruned type
    deleted::Bool               # if the node is supposed to be deleted
    # objs::Vector{JuMP.GenericAffExpr}       # TODO : useless ? erase memory 
    cuts_ref::Vector{CutScore}
    con_cuts::Vector{ConstraintRef}             # TODO : erase memory

    Node() = new()

    """
    A complete node constructor.
    """
    function Node(num::Int64, depth::Int64; pred::Node=Node(), succs::Vector{Node}=Vector{Node}(), var::Int64=0, var_bound::Int64=0)
        n = new()
        n.num = num
        n.depth = depth
    
        n.pred = pred
        n.succs = succs
        n.var = var
        n.var_bound = var_bound
    
        n.RBS = RelaxedBoundSet()
        n.activated = true
        n.pruned = false
        n.prunedType = NONE
        n.deleted = false
        # n.objs = Vector{JuMP.GenericAffExpr}()
        n.cuts_ref = Vector{CutScore}()
        n.con_cuts = Vector{ConstraintRef}()
    
        f(t) = @async println("Finalizing node $(t.num).")
        finalizer(f, n)
    
        # return n
    end
    
end


"""
Return `true` if the given node is the root of a branch-and-bound tree.
"""
function isRoot(node::Node)
    return node.depth == 0 # !isdefined(node, :pred) 
end

"""
Delete the given node in B&B tree. (should be a private function)
"""
function Base.delete!(node::Node)
    node.deleted = true
    node = nothing               # remove from the memory
end

"""
Prune the given node in a B&B tree and delete all successors of the pruned node.
"""
function prune!(node::Node, reason::PrunedType)
    node.pruned = true
    node.prunedType = reason
    to_delete = node.succs[:]
    node.succs = Vector{Node}()

    while length(to_delete) > 0
        n = pop!(to_delete)
        to_delete = vcat(to_delete, n.succs[:])
        delete!(n)
    end
end

"""
From the actual node, go up to the root to get the partial assignment of variables.
"""
function getPartialAssign(actual::Node)
    assignment = Dict{Int64, Int64}() # var index => bound value
    if isRoot(actual) # the actual node is the root 
        return assignment
    end
    predecessor = actual.pred
    assignment[actual.var] = actual.var_bound

    while !isRoot(predecessor)     
        actual = predecessor
        predecessor = actual.pred
        assignment[actual.var] = actual.var_bound
    end
    return assignment
end

"""
Given a patrial assignment of variables, remove the fixed bounding.
"""
function removeBounds(pb::BO01Problem, assignment::Dict{Int64, Int64})
    for (var, bound) in assignment
        if bound == 0
            JuMP.set_upper_bound(pb.varArray[var], 1)
        elseif bound == 1
            JuMP.set_lower_bound(pb.varArray[var], 0)
        end
    end
end

"""
Given a partial assignment on variables values, add the corresponding bounds.
"""
function setBounds(pb::BO01Problem, assignment::Dict{Int64, Int64})
    for (var, bound) in assignment
        if bound == 0
            JuMP.set_upper_bound(pb.varArray[var], 0)
        elseif bound == 1
            JuMP.set_lower_bound(pb.varArray[var], 1)
        end
    end
end


