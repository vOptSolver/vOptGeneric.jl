include("struct.jl")

"""
Definition of the node object in B&B tree.
"""
mutable struct Node
    id::Int64                   # identify number (count)
    depth::Int64                # depth in tree
    pred::Int64                 # predecessor's indexes
    succs::Vector{Int64}        # successors' indices
    var::Int64                  # index of the chosen variable to be split
    var_bound::Int64            # variable bound
    RBS::RelaxedBoundSet        # local relaxed bound set
    activated::Bool               # if the node is active
    isLeaf::Bool               # if the node is leave
    isPruned::Bool              # if the node is pruned
    prunedType::PrunedType      # if the node is fathomed, restore pruned type
end

function Node()
    return Node(0, -1, 0, Vector{Int64}(), 0, 0, RelaxedBoundSet(), true, false, false, NONE)
end

"""
The branch and bound tree.
"""
mutable struct BBTree
    tab::Vector{Node}
end

function BBTree()
    return BBTree(Vector{Node}())
end

"""
Delete the i-th element in B&B tree. (should be a private function)
"""
function Base.delete!(tree::BBTree, i::Int64)
    @assert i ≤ length(tree.tab) "index out of bound !"
    for node in tree.tab[i+1:end] 
        node.id = node.id -1
    end
    for node in tree.tab
        if node.pred > i
            node.pred = node.pred -1
        end
        filter!(c -> c!=i, node.succs) # delete successor i, if the node has
        for j = 1:length(node.succs) 
            if node.succs[j] > i node.succs[j] = node.succs[j]-1 end
        end
    end
    deleteat!(tree.tab, i)
end

"""
Prune the given node in a B&B tree and delete all successors of the pruned node.
"""
function prune!(tree::BBTree, node_id::Int64, reason::PrunedType)
    tree.tab[node_id].isPruned = true
    tree.tab[node_id].prunedType = reason
    to_delete = tree.tab[node_id].succs[:]
    tree.tab[node_id].succs = Vector{Int64}()

    while length(to_delete) > 0
        i = pop!(to_delete)
        to_delete = vcat(to_delete, tree.tab[i].succs[:])
        delete!(tree, i)
    end
end

"""
From the actual node, go up to the root to get the partial assignment of variables.
"""
function getPartialAssign(tree::BBTree, actual::Int64)
    assignment = Dict{Int64, Int64}() # var index => bound value
    if actual == 1 # the actual node is the root 
        return assignment
    end
    predecessor = tree.tab[actual].pred
    assignment[tree.tab[actual].var] = tree.tab[actual].var_bound

    while predecessor > 1       
        actual = predecessor
        predecessor = tree.tab[actual].pred
        assignment[tree.tab[actual].var] = tree.tab[actual].var_bound
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


