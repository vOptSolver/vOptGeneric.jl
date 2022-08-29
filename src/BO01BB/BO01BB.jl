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


"""
Retrieve constraints matrix `A` and right hand side vector `b` in standard LP form `Ax≤b`.
"""
function standard_form(pb::BO01Problem)
    # start = time()

    cstrEqualTo, cstrGreaterThan, cstrLessThan, cstrInterval = [
        JuMP.all_constraints(pb.m, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, set_type) 
            for set_type in (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64})
    ]
    numRows = sum(length, (cstrEqualTo, cstrGreaterThan, cstrLessThan)) ; numRows += 2*length(cstrInterval)
    numVars = length(pb.varArray)
    varIndex = Dict(pb.varArray[i] => i for i=1:length(pb.varArray))

    pb.A = zeros(numRows, numVars)
    pb.b = zeros(numRows)
    pb.c = zeros(2, numVars + 1) # ∑(a_i*x_i) + c

    # objectives 
    vd = getvOptData(pb.m)
    for i=1:2 
        pb.c[i, 1] = vd.objs[i].constant
        for (var, coeff) in vd.objs[i].terms
            pb.c[i, varIndex[var]+1] = coeff
        end
    end

    cstr_index = 0
    
    # Ax=b
    for cstrRef in cstrEqualTo
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.value

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            pb.A[cstr_index, varIndex[term]] = coeff
        end
        pb.b[cstr_index] = rhs
    end

    # Ax≥b
    for cstrRef in cstrGreaterThan
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.lower

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            pb.A[cstr_index, varIndex[term]] = coeff * -1
        end
        pb.b[cstr_index] = rhs * -1
    end

    # Ax≤b
    for cstrRef in cstrLessThan
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.upper

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            pb.A[cstr_index, varIndex[term]] = coeff
        end
        pb.b[cstr_index] = rhs
    end
    

    # lb ≤ AX ≤ ub
    for cstrRef in cstrInterval
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        lb = con.set.lower
        ub = con.set.upper

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            pb.A[cstr_index, varIndex[term]] = coeff ; pb.A[cstr_index+1, varIndex[term]] = coeff * -1
        end
        pb.b[cstr_index] = ub ; pb.b[cstr_index+1] = lb * -1
        cstr_index += 1
    end

    # loadT = round(time() - start, digits = 2)
    # @info "loading matrix ... = $loadT"
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
    # TODO : liberate memory space of the parent node 
    if !isRoot(node)
        if length(node.pred.succs) != 2 || node.pred.succs[1].activated || node.pred.succs[2].activated
            nothing
        else
            node.pred.RBS = RelaxedBoundSet()
            if pb.param.cut_activated
                node.pred.cuts_ref = Vector{CutScore}() ; node.pred.con_cuts = Vector{ConstraintRef}()
                node.pred.cutpool = CutPool()
            end
            # nothing
        end
    end


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
        varArray, m, BBparam(), StatInfo(), Matrix{Float64}(undef, 0,0), Vector{Float64}(), Matrix{Float64}(undef, 0,0)
    )

    if cut
        standard_form(problem)
        problem.param.cut_activated = cut ; problem.info.cuts_activated = cut 
    end
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

    problem.info.cuts_infos.cuts_total = problem.info.cuts_infos.cuts_applied
    println("\n total cuts : ", problem.info.cuts_infos.cuts_applied)

    # for (k,v) in problem.cpool.hashMap
    #     println(k, " => ", size(v, 1))
    # end
    return problem.info
end
