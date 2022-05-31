# This file contains main functions for the bi-objective 0-1 branch and bound algorithm.
include("../vOptData.jl")
include("BBtree.jl")
include("branching.jl")
include("fathoming.jl")


function formatting(m::JuMP.Model)
    converted = false; f = []
    vd = m.ext[:vOpt]::vOptData
    @assert length(vd.objSenses) == 2 "this algorithm requires 2 objective functions !"
    for i = 1:2 
       if vd.objSenses[i] == MOI.MAX_SENSE
            vd.objSenses[i] = MOI.MIN_SENSE
            vd.objs[i] *= -1
            converted = true
            push!(f, i)
       end 
    end
    return converted, f
end

function reversion(m::JuMP.Model, f, incumbent::IncumbentSet)
    vd = m.ext[:vOpt]::vOptData
    for i in f
        vd.objSenses[i] = MOI.MAX_SENSE
        vd.objs[i] *= -1
    end

    for i=1:length(incumbent.natural_order_vect) 
        incumbent.natural_order_vect.sols[i].y *= -1
    end
end

"""
The BO01B&B procedure at each iteration. 
Argument :
    - ind : actual node's index in tree tableau
    - pb : BO01Problem 
"""
function iterative_procedure(todo, node::Node, pb::BO01Problem, incumbent::IncumbentSet, round_results, verbose; args...)
    if verbose
        @info "we are et node $(node.num) !"
    end
    # get the actual node
    @assert node.activated == true "the actual node is not activated "
    node.activated = false

    #--------------------
    # test dominance 
    #--------------------
    start = time()
    if fullyExplicitDominanceTest(node, incumbent)
        prune!(node, DOMINANCE)
        if verbose
            @info "node $(node.num) is fathomed by dominance !"
        end
        pb.info.nb_nodes_pruned += 1
        pb.info.test_dom_time += (time() - start)
        return
    end
    pb.info.test_dom_time += (time() - start)


    #-----------------------------------------
    # branching variable + generate new nodes
    #-----------------------------------------
    var_split = pickUpAFreeVar(node, pb)
    if var_split == 0 return end       # is a leaf

    node1 = Node(
        pb.info.nb_nodes + 1, node.depth + 1, 
        pred = node,
        var = var_split, var_bound = 1
    )
    pb.info.nb_nodes += 1

    if LPRelaxByDicho(node1, pb, round_results, verbose; args...) || updateIncumbent(node1, pb, incumbent, verbose)
        node1.activated = false
    else
        addTodo(todo, pb, node1)
    end

    node2 = Node(
        pb.info.nb_nodes + 1, node.depth + 1,
        pred = node, 
        var = var_split, var_bound = 0
    )
    pb.info.nb_nodes += 1

    if LPRelaxByDicho(node2, pb, round_results, verbose; args...) || updateIncumbent(node2, pb, incumbent, verbose)
        node2.activated = false
    else
        addTodo(todo, pb, node2)
    end

    node.succs = [node1, node2]
end

function post_processing(m::JuMP.Model, problem::BO01Problem, incumbent::IncumbentSet, round_results, verbose; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)

    for sol in incumbent.natural_order_vect.sols
        for x in sol.xEquiv
            push!(vd.Y_N, round_results ? round.(sol.y) : sol.y)
            push!(vd.X_E, x)
        end
    end

    problem.info.relaxation_time = round(problem.info.relaxation_time, digits = 2)
    problem.info.test_dom_time = round(problem.info.test_dom_time, digits = 2)
    problem.info.update_incumb_time = round(problem.info.update_incumb_time, digits = 2)
end

"""
A bi-objective binary(0-1) branch and bound algorithm.
"""
function solve_branchbound(m::JuMP.Model, round_results, verbose; args...)
    start = time()

    converted, f = formatting(m)

    varArray = JuMP.all_variables(m)
    problem = BO01Problem(varArray, m, BBparam(), StatInfo())

    # initialize the incumbent list by heuristics or with Inf
    incumbent = IncumbentSet()

    # by default, we take the breadth-first strategy (FIFO queue)
    todo = initQueue(problem)

    # step 0 : create the root and add to the todo list
    root = Node(problem.info.nb_nodes +1, 0)
    problem.info.nb_nodes += 1

    if LPRelaxByDicho(root, problem, round_results, verbose; args...) || updateIncumbent(root, problem, incumbent, verbose)
        if converted
            reversion(m, f, incumbent)
        end
        post_processing(m, problem, incumbent, round_results, verbose; args...)
        problem.info.total_times = round(time() - start, digits = 2)
        return problem.info
    end

    addTodo(todo, problem, root)

    # continue to fathom the node until todo list is empty
    while length(todo) > 0
        node_ref = nextTodo(todo, problem)
        
        if node_ref[].deleted
            finalize(node_ref[])
        else
            iterative_procedure(todo, node_ref[], problem, incumbent, round_results, verbose; args...)
        end
    end

    if converted
        reversion(m, f, incumbent)
    end
    post_processing(m, problem, incumbent, round_results, verbose; args...)
    problem.info.total_times = round(time() - start, digits = 2)
    return problem.info
end