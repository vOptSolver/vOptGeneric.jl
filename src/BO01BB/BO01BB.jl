# This file contains main functions for the bi-objective 0-1 branch and bound algorithm.

include("struct.jl")
include("../algorithms.jl")

"""
The BO01B&B procedure at each iteration.

Argument :
    - node : actual node defining a subproblem
    - pb : BO01Problem 
"""
function iterative_procedure(ind::Int64, BB_tree::BBTree, pb::BO01Problem, incumbent::Incumbent, round_results, verbose; args...)
    # get the actual node
    node = BB_tree.tree[ind]
    node.actived = false
    if node.var > 0
        node.assignment[node.var] = node.var_bound
    end
    
    #-----------------------------------------------------------
    # STEP 1 : resolving the LP relaxation + check feasiblility
    #-----------------------------------------------------------
    undo_relax = relax_integrality(pb.m)
    setBounds(pb, node.assignment)
    solve_dicho(pb.m, round_results, verbose ; args...)
    removeBounds(pb, node.assignment)
    undo_relax()
    vd_LP = getvOptData(pb.m)

    if size(vd_LP.Y_N, 1) == 0
        to_release = prune!(node, INFEASIBILITY)
        release(BB_tree, to_release)
        return 
    end

    # construct/complete the relaxed bound set
    for i = 1:length(vd_LP.Y_N)
        push!(node.RBS.sols, Sol(vd_LP.X_E[i], vd_LP.Y_N[i]))
    end

    #-----------------------------------------------------------
    # STEP 2 : check optimality && update the incumbent set
    #-----------------------------------------------------------
    all_integer = true
    for i = 1:length(node.RBS.sols)
        if node.RBS.sols.sol_vect[i].is_integer
            s = node.RBS.sols.sol_vect[i]
            push!(incumbent.sols, s)
        else
            all_integer = false
        end
    end

    if all_integer
        to_release = prune!(node, OPTIMALITY)
        release(BB_tree, to_release)
        return
    end

    #---------------------------
    # STEP 3 : test dominance 
    #---------------------------

    #------------------------------
    # STEP 4 : branching variable 
    #------------------------------
    free_vars = [ind for ind in 1:length(pb.varArray) if node.assignment[ind] == -1]
    var_split = free_vars[rand(1:length(free_vars))]
    node1 = Node(
        length(BB_tree.tree) + 1,
        node.id,
        Vector{Int64}(),
        var_split, 1,
        RelaxedBoundSet(), NatrualOrderVector(),
        true, false, NONE, node.assignment[:]
    )
    push!(BB_tree.tree, node1)

    node2 = Node(
        length(BB_tree.tree) + 1,
        node.id,
        Vector{Int64}(),
        var_split, 0,
        RelaxedBoundSet(), NatrualOrderVector(),
        true, false, NONE, node.assignment[:]
    )
    push!(BB_tree.tree, node2)

    node.succs = [node1.id, node2.id]

    return 
end



"""
A bi-objective binary(0-1) branch and bound algorithm.
"""
function solve_branchbound(m::JuMP.Model, round_results, verbose; args...)

    println("hello MO01BB :)")

    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    varArray = JuMP.all_variables(m)
    problem = BO01Problem(vd, varArray, m, BBparam())

    # initialize the incumbent list by heuristics or with Inf
    incumbent = IncumbentSet()

    BB_tree = BBTree()

    # by defaut, we take the breadth-first strategy (FIFO queue)
    todo = Queue{Int64}()

    # step 0 : create the root and add to the todo list
    root = Node(
        length(BB_tree.tree) + 1,
        0, 
        Vector{Int64}(),
        0, 0,
        RelaxedBoundSet(), NatrualOrderVector(),
        true, false, NONE, 
        [-1 for _ in 1:length(varArray)]
    )
    push!(BB_tree.tree, root)
    enqueue!(todo, root.id)

    # continue to fathom the node until todo list is empty
    while length(todo) > 0
        ind = dequeue!(todo)
        iterative_procedure(ind, BB_tree, problem, incumbent, round_results, verbose; args...)

        # add successors of node at tree[ind] to the waitting queue
        for id in BB_tree.tree[ind].succs
            enqueue!(todo, id)
        end

    end

end