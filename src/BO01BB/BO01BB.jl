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
function iterative_procedure(todo, ind::Int64, tree::BBTree, pb::BO01Problem, incumbent::IncumbentSet, round_results, verbose; args...)
    if verbose
        @info "we are at node ", ind
    end
    # get the actual node
    node = tree.tab[ind]
    @assert node.activated == true "the actual node is not activated "
    node.activated = false

    #---------------------------
    # STEP 3 : test dominance 
    #---------------------------

    #-----------------------------------------
    # branching variable + generate new nodes
    #-----------------------------------------
    var_split = pickUpAFreeVar(ind, tree, pb)
    if tree.tab[ind].isLeaf
        return
    end

    node1 = Node(
        length(tree.tab) + 1,
        tree.tab[ind].depth + 1, ind,
        Vector{Int64}(),
        var_split, 1,
        RelaxedBoundSet(),
        true, false, false, NONE
    )
    push!(tree.tab, node1)
    if LPRelaxByDicho(node1.id, tree, pb, round_results, verbose; args...) || updateIncumbent(node1.id, tree, incumbent, verbose)
        node1.activated = false
    else
        addTodo(todo, pb, node1.id)
    end

    node2 = Node(
        length(tree.tab) + 1,
        tree.tab[ind].depth + 1, ind,
        Vector{Int64}(),
        var_split, 0,
        RelaxedBoundSet(),
        true, false, false, NONE
    )
    push!(tree.tab, node2)
    if LPRelaxByDicho(node2.id, tree, pb, round_results, verbose; args...) || updateIncumbent(node2.id, tree, incumbent, verbose)
        node2.activated = false
    else
        addTodo(todo, pb, node2.id)
    end

    tree.tab[ind].succs = [node1.id, node2.id]
end



"""
A bi-objective binary(0-1) branch and bound algorithm.
"""
function solve_branchbound(m::JuMP.Model, round_results, verbose; args...)

    println("hello BO01BB :)")

    converted, f = formatting(m)

    varArray = JuMP.all_variables(m)
    problem = BO01Problem(varArray, m, BBparam(), StatInfo())

    # initialize the incumbent list by heuristics or with Inf
    incumbent = IncumbentSet()

    tree = BBTree()

    # by default, we take the breadth-first strategy (FIFO queue)
    todo = initQueue(problem)

    # step 0 : create the root and add to the todo list
    root = Node(
        1, 0, 0, 
        Vector{Int64}(),
        0, 0,
        RelaxedBoundSet(),
        true, false, false, NONE
    )
    push!(tree.tab, root)

    if LPRelaxByDicho(root.id, tree, problem, round_results, verbose; args...) || updateIncumbent(root.id, tree, incumbent, verbose)
        return
    end

    addTodo(todo, problem, root.id)

    # continue to fathom the node until todo list is empty
    while length(todo) > 0
        println("new iteration , ", todo)
        ind = nextTodo(todo, problem)
        iterative_procedure(todo, ind, tree, problem, incumbent, round_results, verbose; args...)
    end

    if converted
        reversion(m, f, incumbent)
    end

    println("incumbent : ", incumbent.natural_order_vect)
end