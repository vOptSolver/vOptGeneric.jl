# This file contains main functions for the bi-objective 0-1 branch and bound algorithm.
include("../vOptData.jl")
include("BBtree.jl")
include("branching.jl")
include("fathoming.jl")
include("displayGraphic.jl")

using TimerOutputs
tmr = TimerOutput()


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

# TODO : retrieve matrices
"""
Retrieve constraints matrix `A` and right hand side vector `b` in standard LP form `Ax≤b`.
"""
function standard_form(pb::BO01Problem)

    cstrEqualTo, cstrGreaterThan, cstrLessThan, cstrInterval = [
        JuMP.all_constraints(pb.m, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, set_type) 
            for set_type in (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64})
    ]
    numRows = sum(length, (cstrEqualTo, cstrGreaterThan, cstrLessThan))
    numRows += 2*length(cstrInterval)
    println(pb.varArray)
    println("numRows = ", numRows)

    numVars = length(pb.varArray)
    println("numVars = ", numVars)

    A = zeros(numRows, numVars)
    b = zeros(numRows)
    c = zeros(2, numVars + 1) # ∑(a_i*x_i) + c

    # objectives 
    vd = getvOptData(pb.m)
    for i=1:2 
        c[i, 1] = vd.objs[i].constant
        for (var, coeff) in vd.objs[i].terms
            for j=1:numVars
                if var == pb.varArray[j] 
                    c[i, j+1] = coeff ; break
                end 
            end
        end
    end

    cstr_index = 0
    
    # Ax=b
    println("Ax=b")
    for cstrRef in cstrEqualTo
        println(cstrRef)
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.value

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            for i=1:numVars 
                if term == pb.varArray[i]
                    A[cstr_index, i] = coeff
                end
            end
        end
        b[cstr_index] = rhs
    end

    # Ax≥b
    println("Ax≥b")
    for cstrRef in cstrGreaterThan
        println(cstrRef)
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.lower

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            for i=1:numVars 
                if term == pb.varArray[i]
                    A[cstr_index, i] = coeff * -1
                end
            end
        end
        b[cstr_index] = rhs * -1
    end

    # Ax≤b
    println("Ax≤b")
    for cstrRef in cstrLessThan
        println(cstrRef)
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.upper

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            for i=1:numVars 
                if term == pb.varArray[i]
                    A[cstr_index, i] = coeff
                end
            end
        end
        b[cstr_index] = rhs
    end
    

    # lb ≤ AX ≤ ub
    println("lb ≤ AX ≤ ub")
    for cstrRef in cstrInterval
        println(cstrRef)
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        lb = con.set.lower
        ub = con.set.upper

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            for i=1:numVars 
                if term == pb.varArray[i]
                    A[cstr_index, i] = coeff
                    A[cstr_index+1, i] = coeff * -1
                end
            end
        end
        b[cstr_index] = ub
        b[cstr_index+1] = lb * -1
        cstr_index += 1
    end

    println("A : ", A)
    println("b : ", b)
    println("c : ", c)

    pb.A = deepcopy(A)
    pb.b = deepcopy(b)
    pb.c = deepcopy(c)
end

"""
The BO01B&B procedure at each iteration. 
Argument :
    - ind : actual node's index in tree tableau
    - pb : BO01Problem 
"""
function iterative_procedure(todo, node::Node, pb::BO01Problem, incumbent::IncumbentSet, round_results, verbose; args...)
    if verbose
        @info "at node $(node.num) |Y_N| = $(length(incumbent.natural_order_vect))"
    end
    # get the actual node
    @assert node.activated == true "the actual node is not activated "
    node.activated = false

    #--------------------
    # test dominance 
    #--------------------
    start = time()
    if ( @timeit tmr "dominance" fullyExplicitDominanceTest(node, incumbent) )
        prune!(node, DOMINANCE)
        if verbose
            @info "node $(node.num) is fathomed by dominance ! |LBS|=$(length(node.RBS.natural_order_vect))" 
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

    if ( @timeit tmr "relax" LPRelaxByDicho(node1, pb, round_results, verbose; args...) ) || 
        ( @timeit tmr "incumbent" updateIncumbent(node1, pb, incumbent, verbose) )
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

    if ( @timeit tmr "relax" LPRelaxByDicho(node2, pb, round_results, verbose; args...) ) || 
        ( @timeit tmr "incumbent" updateIncumbent(node2, pb, incumbent, verbose) )
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
        push!(vd.Y_N, round_results ? round.(sol.y) : sol.y)
        for x in sol.xEquiv
            push!(vd.X_E, x)
        end
    end

    # for sol in incumbent.natural_order_vect.sols
    #     push!(vd.Y_N, round_results ? round.(sol.y) : sol.y)
    #     push!(vd.X_E, sol.xEquiv[1])
    # end

    s = sortperm(vd.Y_N, by = first)
    vd.X_E, vd.Y_N = vd.X_E[s], vd.Y_N[s]

    problem.info.relaxation_time = round(problem.info.relaxation_time, digits = 2)
    problem.info.test_dom_time = round(problem.info.test_dom_time, digits = 2)
    problem.info.update_incumb_time = round(problem.info.update_incumb_time, digits = 2)

    if problem.info.cuts_activated
        problem.info.cuts_infos.times_calling_dicho = round(problem.info.cuts_infos.times_calling_dicho, digits = 2)
        problem.info.cuts_infos.times_calling_separators = round(problem.info.cuts_infos.times_calling_separators, digits = 2)
        problem.info.cuts_infos.times_oper_cutPool = round(problem.info.cuts_infos.times_oper_cutPool, digits = 2)
        problem.info.cuts_infos.times_total_for_cuts = round(problem.info.cuts_infos.times_total_for_cuts, digits = 2)
        problem.info.cuts_infos.times_add_retrieve_cuts = round(problem.info.cuts_infos.times_add_retrieve_cuts, digits = 2)
    end
end

"""
A bi-objective binary(0-1) branch and bound algorithm.
"""
function solve_branchboundcut(m::JuMP.Model, cut::Bool, round_results, verbose; args...)
    start = time()

    converted, f = formatting(m)

    varArray = JuMP.all_variables(m)
    problem = BO01Problem(
        varArray, m, BBparam(), StatInfo(), Matrix{Float64}(undef, 0,0), Vector{Float64}(), Matrix{Float64}(undef, 0,0), CutPool()
    )

    standard_form(problem)
    problem.param.cut_activated = cut ; problem.info.cuts_activated = cut 

    # relaxation LP
    undo_relax = JuMP.relax_integrality(problem.m)

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
        undo_relax()
        return problem.info
    end

    addTodo(todo, problem, root)

    # continue to fathom the node until todo list is empty
    @timeit tmr "BB loop" begin
        while length(todo) > 0
            @timeit tmr "next node" node_ref = nextTodo(todo, problem)

            if node_ref[].deleted
                finalize(node_ref[])
            else
                @timeit tmr "iteration" iterative_procedure(todo, node_ref[], problem, incumbent, round_results, verbose; args...)
            end
        end
    end

    if converted
        reversion(m, f, incumbent)
    end
    post_processing(m, problem, incumbent, round_results, verbose; args...)
    problem.info.total_times = round(time() - start, digits = 2)
    
    MB = 10^6

    problem.info.tree_size = round(Base.summarysize(root)/MB, digits = 3)
    
    undo_relax()
    show(tmr)

    problem.info.cuts_infos.cuts_total = length(problem.cpool.tab)
    println("\n total cuts : ", length(problem.cpool.tab))
    # for cut in problem.cpool.tab
    #     println("cut : ", cut)
    # end
    return problem.info
end
