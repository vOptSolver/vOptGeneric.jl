# This file contains functions related to node fathoming.

include("BBtree.jl")
include("../algorithms.jl")

"""
Compute and stock the relaxed bound set (i.e. the LP relaxation) of the (sub)problem defined by the given node.

Return `true` if the node is pruned by infeasibility.
"""
function LPRelaxByDicho(node_id::Int64, tree::BBTree, pb::BO01Problem, round_results, verbose ; args...)
    #------------------------------------------------------------------------------
    # solve the LP relaxation by dichotomy method including the partial assignment
    #------------------------------------------------------------------------------
    undo_relax = JuMP.relax_integrality(pb.m)
    assignment = getPartialAssign(tree, node_id) 
    if verbose
        println("LP relax : at node ", node_id, " assignment is ", assignment)
    end
    setBounds(pb, assignment)
    solve_dicho(pb.m, round_results, verbose ; args...)
    removeBounds(pb, assignment)
    undo_relax()
    vd_LP = getvOptData(pb.m)

    #-------------------------------------------------------------------------------
    # in case of the LP relaxed (sub) problem is infeasible, prune the actual node
    #-------------------------------------------------------------------------------
    if size(vd_LP.Y_N, 1) == 0
        prune!(tree, node_id, INFEASIBILITY)
        if verbose
            @info "node ", node_id, "is unfeasible !"
        end
        return true
    end

    # construct/complete the relaxed bound set
    for i = 1:length(vd_LP.Y_N)
        s = Solution(vd_LP.X_E[i], vd_LP.Y_N[i])
        push!(tree.tab[node_id].RBS.natural_order_vect, s)
        if i < length(vd_LP.Y_N)
            tree.tab[node_id].RBS.segments[s] = true
        end
    end 
    if verbose
        println("RBS is ", tree.tab[node_id].RBS.natural_order_vect)
    end
    return false
end


"""
At the given node, update (filtered by dominance) the global incumbent set. 

Return `true` if the node is pruned by optimality.
"""
function updateIncumbent(node_id::Int64, tree::BBTree, incumbent::IncumbentSet, verbose)
    #-----------------------------------------------------------
    # check optimality && update the incumbent set
    #-----------------------------------------------------------
    node = tree.tab[node_id]
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
        prune!(tree, node_id, OPTIMALITY)
        if verbose
            @info "node ", node_id, "is fathomed by optimality !"
        end
        return true
    end
    return false
end


"""
A fully explicit dominance test, and prune the given node if it's fathomed by dominance.

"""
function fullyExplicitDominanceTest(node_id::Int64, tree::BBTree, incumbent::IncumbentSet)
    @assert length(tree.tab[node_id].RBS.natural_order_vect) > 0 "relaxed bound set is empty for node $node_id"

    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(incumbent.natural_order_vect) == 0 return end

    # when the LBS consists of a single point 

    # 
end