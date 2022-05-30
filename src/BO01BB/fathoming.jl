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
    undo_relax = JuMP.relax_integrality(pb.m)
    assignment = getPartialAssign(node) 
    setBounds(pb, assignment)
    solve_dicho(pb.m, round_results, false ; args...)
    removeBounds(pb, assignment)
    undo_relax()
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
        pb.info.relaxation_time += round(time() - start, digits = 2)
        return true
    end

    # construct/complete the relaxed bound set
    for i = 1:length(vd_LP.Y_N)
        push!(node.RBS.natural_order_vect, Solution(vd_LP.X_E[i], vd_LP.Y_N[i]))
    end 
    for i=1:length(node.RBS.natural_order_vect)-1
        node.RBS.segments[node.RBS.natural_order_vect.sols[i]] = true
    end

    pb.info.relaxation_time += round(time() - start, digits = 2)
    return false
end


"""
At the given node, update (filtered by dominance) the global incumbent set. 

Return `true` if the node is pruned by optimality.
"""
function updateIncumbent(node::Node, incumbent::IncumbentSet, verbose)
    start = time()
    #-----------------------------------------------------------
    # check optimality && update the incumbent set
    #-----------------------------------------------------------
    all_binary = true

    for i = 1:length(node.RBS.natural_order_vect)
        if node.RBS.natural_order_vect.sols[i].is_binary
            s = node.RBS.natural_order_vect.sols[i]
            push!(incumbent.natural_order_vect, s, filtered=true)
        else
            all_binary = false
        end
    end

    if all_binary
        prune!(node, OPTIMALITY)
        if verbose
            @info "node $(node.num) is fathomed by optimality !"
        end
        pb.info.nb_nodes_pruned += 1
        pb.info.update_incumb_time += round(time() - start, digits = 2)
        return true
    end
    pb.info.update_incumb_time += round(time() - start, digits = 2)
    return false
end

"""
Return a list of local nadir points (so-called cornor points) of the given incumbent set, or the single point it contains.
"""
function getNadirPoints(incumbent::IncumbentSet)
    if length(incumbent.natural_order_vect) == 1
        return [incumbent.natural_order_vect.sols[1]]
    end
    nadir_pts = Vector{Solution}()
    for i = 1:length(incumbent.natural_order_vect)-1 
        push!(nadir_pts, Solution(
            Vector{Vector{Float64}}(),
            [max(incumbent.natural_order_vect.sols[i].y[1], incumbent.natural_order_vect.sols[i+1].y[1]), 
            max(incumbent.natural_order_vect.sols[i].y[2], incumbent.natural_order_vect.sols[i+1].y[2])],
            false )
        )
    end

    return nadir_pts
end

"""
A fully explicit dominance test, and prune the given node if it's fathomed by dominance.
(i.e. ∀ l∈L: ∃ u∈U s.t. λu ≤ λl )

Return `true` if the given node is fathomed by dominance.
"""
function fullyExplicitDominanceTest(node::Node, incumbent::IncumbentSet)
    @assert length(node.RBS.natural_order_vect) > 0 "relaxed bound set is empty for node $node_id"

    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(incumbent.natural_order_vect) == 0 return false end

    # if there exists an upper bound u s.t. u≦l
    function weak_dom(l)
        for u in incumbent.natural_order_vect.sols
            if u <= l
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
    # if the LBS consists of segments and points
    # ----------------------------------------------
    nadir_pts = getNadirPoints(incumbent)

    for i=1:length(node.RBS.natural_order_vect)-1
        sol = node.RBS.natural_order_vect.sols[i]

        # in case of segment
        if node.RBS.segments[sol]
            sol2 = node.RBS.natural_order_vect.sols[i+1]
            λ = [abs(sol.y[2] - sol2.y[2]), abs(sol.y[1] - sol2.y[1])]      # normal to the segment

            for u in incumbent.natural_order_vect.sols
                if λ'*u.y <= λ'*sol2.y
                    return true
                end
            end

        # in case of isolated point
        elseif weak_dom(sol)
            return true
        end
    end

    return false
end