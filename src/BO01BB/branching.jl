# This file contains strategies on variables splitting and tree traversing.

using DataStructures

include("struct.jl")

"""
Return an initialized todo list according to the fixed parameter.
"""
function initQueue(pb::BO01Problem)
    if pb.param.traverse == :bfs
        return Queue{Int64}()
    elseif pb.param.traverse == :dfs
        return Stack{Int64}() 
    elseif pb.param.traverse == :arbitrary
        return Vector{Int64}()
    else
        @error "Unknown traverse parameter !"
    end
end


"""
Add a node identify in the todo list.
"""
function addTodo(todo, pb::BO01Problem, node_id::Int64)
    if pb.param.traverse == :bfs
        enqueue!(todo, node_id)
    elseif pb.param.traverse == :dfs
        push!(todo, node_id) 
    elseif pb.param.traverse == :arbitrary
        push!(todo, node_id)
    else
        @error "Unknown traverse parameter !"
    end
end

"""
Return the next element in the todo list.
"""
function nextTodo(todo, pb::BO01Problem)
    if pb.param.traverse == :bfs
        return dequeue!(todo)
    elseif pb.param.traverse == :dfs
        return pop!(todo) 
    elseif pb.param.traverse == :arbitrary
        i = rand(1:length(todo))
        next = todo[i]
        deleteat!(todo, i)
        return next
    else
        @error "Unknown traverse parameter !"
    end
end

"""
Pich up a free variable to be split according to the prefiexd strategy.
"""
function pickUpAFreeVar(actual::Int64, tree::BBTree, pb::BO01Problem)
    if pb.param.branching == :arbitrary
        free_vars = [ind for ind in 1:length(pb.varArray)]
        fixed_var = collect(keys(getPartialAssign(tree, actual)))
        filter!(v -> v ∉ fixed_var, free_vars)
        if length(free_vars) == 0
            tree.tab[actual].isLeaf = true
            return 0
        else
            return free_vars[rand(1:length(free_vars))]
        end
    else
        @error "Unknown branching parameter !"
    end
end