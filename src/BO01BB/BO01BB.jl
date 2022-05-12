# This file contains main functions for the bi-objective 0-1 branch and bound algorithm.

include("struct.jl")
include("../algorithms.jl")

"""
The BO01B&B procedure at each iteration.

Argument :
    - node : actual node defining a subproblem
    - m : JuMP model 

Return : 
    - node1 : the new generated subproblem with a variable bounded
    - node2 : the new generated subproblem with a variable bounded complementary
"""
function iterative_procedure(node::Node, m::JuMP.Model, round_results, verbose; args...)
    node.actived = false

    # STEP 1 : resolving the LP relaxation + check feasiblility
    undo_relax = relax_integrality(m)
    solve_dicho(m, round_results, verbose ; args...)
    undo_relax()
    vd = getvOptData(m)

    if size(vd.Y_N, 1) == 0
        prune!(node, INFEASIBILITY)
        return
    end

    # STEP 2 : check optimality

    # STEP 3 : test dominance 

    # STEP 4 : branching variable 

    return node1, node2
end



"""
A bi-objective binary(0-1) branch and bound algorithm.
"""
function solve_branchbound(m::JuMP.Model, round_results, verbose; args...)

    println("hello MO01BB :)")

    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    varArray = JuMP.all_variables(m)
    problem = BO01Problem(vd, varArray, m)

    # initialize the incumbent list by heuristics or with Inf
    incumbent = IncumbentSet()

    BB_tree = Vector{Node}()

    # by defaut, we take the breadth-first strategy (FIFO queue)
    todo = Queue{Int64}()

    # step 0 : create the root and add to the todo list
    root = Node()
    push!(BB_tree, root)
    enqueue!(todo, root.id)

    # step 1 : continue to fathom the node until todo list is empty
    while length(todo) > 0
        ind = dequeue!(todo)
        node1, node2 = iterative_procedure(BB_tree[ind], m, round_results, verbose; args...)
        push!(BB_tree, node1); push!(BB_tree, node2)
        enqueue!(todo, node1.id); enqueue!(todo, node2.id)
    end

end