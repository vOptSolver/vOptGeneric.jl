# This file contains strategies on variables splitting and tree traversing.

using DataStructures

include("BBtree.jl")

"""
Return an initialized todo list according to the fixed parameter.
"""
function initQueue(pb::BO01Problem)
    if pb.param.traverse == :bfs
        return Queue{Base.RefValue{Node}}()
    elseif pb.param.traverse == :dfs
        return Stack{Base.RefValue{Node}}() 
    elseif pb.param.traverse == :arbitrary
        return Vector{Base.RefValue{Node}}()
    else
        @error "Unknown traverse parameter !"
    end
end


"""
Add a node identify in the todo list.
"""
function addTodo(todo, pb::BO01Problem, node::Node)
    if pb.param.traverse == :bfs
        enqueue!(todo, Ref(node))
    elseif pb.param.traverse == :dfs
        push!(todo, Ref(node)) 
    elseif pb.param.traverse == :arbitrary
        push!(todo, Ref(node))
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
function pickUpAFreeVar(assignment::Dict{Int64, Int64}, pb::BO01Problem)
    if pb.param.branching == :arbitrary
        free_vars = [ind for ind in 1:length(pb.varArray)]
        fixed_var = collect(keys(assignment))
        filter!(v -> v ∉ fixed_var, free_vars)
        return (length(free_vars) > 0) ? free_vars[rand(1:length(free_vars))] : 0
    else
        @error "Unknown branching parameter !"
    end
end