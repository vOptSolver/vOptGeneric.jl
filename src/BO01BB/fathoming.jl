# This file contains functions related to node fathoming.

include("BBtree.jl")
include("../algorithms.jl")

"""
Compute and stock the relaxed bound set (i.e. the LP relaxation) of the (sub)problem defined by the given node.
Return `true` if the node is pruned by infeasibility.
"""
function LPRelaxByDicho(node::Node, pb::BO01Problem, round_results, verbose ; args...)
    start = time()
    #------------------------------------------------------------------------------
    # solve the LP relaxation by dichotomy method including the partial assignment
    #------------------------------------------------------------------------------
    # undo_relax = JuMP.relax_integrality(pb.m)
    assignment = getPartialAssign(node)
    setBounds(pb, assignment)
    solve_dicho(pb.m, round_results, false ; args...)
    removeBounds(pb, assignment)
    # undo_relax()

    vd_LP = getvOptData(pb.m)

    #-------------------------------------------------------------------------------
    # in case of the LP relaxed (sub) problem is infeasible, prune the actual node
    #-------------------------------------------------------------------------------
    if size(vd_LP.Y_N, 1) == 0
        prune!(node, INFEASIBILITY)
        if verbose
            @info "node $(node.num) is unfeasible !"
        end
        pb.info.nb_nodes_pruned += 1
        pb.info.relaxation_time += (time() - start)
        # pb.info.status = MOI.INFEASIBLE
        return true
    end

    # construct/complete the relaxed bound set
    for i = 1:length(vd_LP.Y_N)
        push!(node.RBS.natural_order_vect, Solution(vd_LP.X_E[i], vd_LP.Y_N[i]))
    end
    for i=1:length(node.RBS.natural_order_vect)-1
        node.RBS.segments[node.RBS.natural_order_vect.sols[i]] = true
    end

    pb.info.relaxation_time += (time() - start)
    return false
end


"""
At the given node, update (filtered by dominance) the global incumbent set.
Return `true` if the node is pruned by optimality.
"""
function updateIncumbent(node::Node, pb::BO01Problem, incumbent::IncumbentSet, verbose)
    start = time()
    #-----------------------------------------------------------
    # check optimality && update the incumbent set
    #-----------------------------------------------------------

    for i = 1:length(node.RBS.natural_order_vect)
        if node.RBS.natural_order_vect.sols[i].is_binary
            s = node.RBS.natural_order_vect.sols[i]
            push!(incumbent.natural_order_vect, s, filtered=true)
        end
    end

    if length(node.RBS.natural_order_vect)==1 && node.RBS.natural_order_vect.sols[1].is_binary
        prune!(node, OPTIMALITY)
        if verbose
            @info "node $(node.num) is fathomed by optimality ! and length = $(length(node.RBS.natural_order_vect))"
        end
        pb.info.nb_nodes_pruned += 1
        pb.info.update_incumb_time += (time() - start)
        return true
    end
    pb.info.update_incumb_time += (time() - start)
    return false
end

"""
Return local nadir points (so-called corner points) of the given incumbent set, or the single point it contains.
"""
function getNadirPoints(incumbent::IncumbentSet, ptl, ptr)
    nadir_pts = NaturalOrderVector()
    @assert length(incumbent.natural_order_vect) > 1 "`getNadirPoints` requires at least two upper bounds in incumbent list."

    # condition ideal point
    if incumbent.natural_order_vect.sols[1].y[1] <= ptr.y[1] && incumbent.natural_order_vect.sols[1].y[2] <= ptl.y[2]
        return (nadir_pts, true)
    end

    for i = 1:length(incumbent.natural_order_vect)-1
        # condition ideal point
        if incumbent.natural_order_vect.sols[i+1].y[1] <= ptr.y[1] && incumbent.natural_order_vect.sols[i+1].y[2] <= ptl.y[2]
            return (nadir_pts, true)
        end

        push!(nadir_pts, Solution(
            Vector{Vector{Float64}}(),
            [incumbent.natural_order_vect.sols[i].y[1],
            incumbent.natural_order_vect.sols[i+1].y[2]
            ],
            true ) , filtered=true
        )
    end

    return (nadir_pts, false)
end

"""
A fully explicit dominance test, and prune the given node if it's fathomed by dominance.
(i.e. ∀ l∈L: ∃ u∈U s.t. λu ≤ λl )
Return `true` if the given node is fathomed by dominance.
"""
function fullyExplicitDominanceTest(node::Node, incumbent::IncumbentSet)
    @assert length(node.RBS.natural_order_vect) > 0 "relaxed bound set is empty for node $(node.num)"

    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(incumbent.natural_order_vect) == 0 return false end

    # if there exists an upper bound u s.t. u≦l
    function weak_dom(l)
        for u ∈ incumbent.natural_order_vect.sols
            if u ≤ l
                return true
            end
        end
        return false
    end

    # ------------------------------------------
    # if the LBS consists of a single point
    # ------------------------------------------
    if length(node.RBS.natural_order_vect) == 1
        l = node.RBS.natural_order_vect.sols[1]
        return weak_dom(l)
    end

    # ----------------------------------------------
    # if the LBS consists of segments
    # ----------------------------------------------
    if length(incumbent.natural_order_vect) == 1 return false end

    #ideal point of LBS
    ptl = node.RBS.natural_order_vect.sols[1]
    ptr = node.RBS.natural_order_vect.sols[end]

    (nadir_pts, fathomed) = getNadirPoints(incumbent, ptl, ptr)
    if fathomed return true end

    # test range condition necessary
    # if length(nadir_pts) == 1 return fathomed end

    u_l = incumbent.natural_order_vect.sols[1]
    u_r = incumbent.natural_order_vect.sols[end]

    sufficient = (u_l.y[2] < ptl.y[2] && u_r.y[1] < ptr.y[1])

    if !sufficient return false end

    # iterate of all local nadir points
    for u ∈ nadir_pts.sols
        existence = false
        compared = false

        # condition ideal point
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        end

        # a complete pairwise comparison
        for i=1:length(node.RBS.natural_order_vect)-1              # ∀ segment l ∈ LBS 

            sol_l = node.RBS.natural_order_vect.sols[i]
            sol_r = node.RBS.natural_order_vect.sols[i+1]

            if (u.y[1] > sol_l.y[1] || u.y[1] < sol_r.y[1]) && (u.y[2] > sol_r.y[2] || u.y[2] < sol_l.y[2])
                continue
            end
            

            λ = [sol_r.y[2] - sol_l.y[2], sol_l.y[1] - sol_r.y[1]]      # normal to the segment

            compared = true

            if λ'*u.y < λ'*sol_r.y #&& λ'*u.y < λ'*sol_l.y
                existence = true
                break
            end
        end
        
        # condition dominance violated
        if compared && !existence return false end

        if !compared && (u.y[1] ≥ ptr.y[1] && u.y[2] ≥ ptl.y[2])
            return false
        end
    end

    return true
end