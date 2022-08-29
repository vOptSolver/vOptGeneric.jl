# This file contains functions related to node fathoming.

include("cuttingPlanes.jl")


"""
Given a point `x_star`, iterate all valid cuts of parent node and stock the references 
in the current node.
"""
function loadingCutInPool(node::Node, pb::BO01Problem)
    if isRoot(node) return end 
    # --------------------------------------------------------------------------
    # iterate in the global cut pool, identify + stock the violated cuts indexes
    # --------------------------------------------------------------------------
    # @info "node $(node.num) loadingCutInPool ... pred = $(node.pred.num)"
    # println("pred.cutpool : ")
    # for (k,v) in node.pred.cutpool.hashMap
    #     println(k, " => ", size(v, 1))
    # end
    # println("pred.cutrefs : \n", node.pred.cuts_ref)


    l = 1 ; LBS = node.RBS.natural_order_vect.sols

    while l ≤ length(LBS)
        if LBS[l].is_binary 
            l += 1 ; continue
        end

        xₗ_star = LBS[l].xEquiv[1]

        for ∇ = max_step:-1:0 
            if ∇ == 0
                # single-point cut 
                for (k, cuts) in node.pred.cutpool.hashMap
                    for α in cuts
                        violationₗ = maximum([ (xₗ_star'*α[2:end] - α[1]), 0.0 ])
                        if violationₗ > 0.0
                            ineq = Cut(α)
                            if push!(node.cutpool, ineq)# && push_cutScore(node.cuts_ref, CutScore(length(node.cutpool.hashMap[k]), violationₗ, k))
                                con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                            end
                        end
                    end
                end
                l += 1
            else 
                applied = false
                r = l+∇
                if r > length(LBS) || LBS[r].is_binary continue end

                xᵣ_star = LBS[r].xEquiv[1]
                # multi-point cut 
                for (k, cuts) in node.pred.cutpool.hashMap
                    for α in cuts 
                        violationₗ = maximum([ (xₗ_star'*α[2:end] - α[1]), 0.0 ])
                        violationᵣ = maximum([ (xᵣ_star'*α[2:end] - α[1]), 0.0 ])
                        viol = maximum([violationₗ, violationᵣ])
                        if viol > 0.0
                            applied = true
                            ineq = Cut(α)
                            if push!(node.cutpool, ineq)# && push_cutScore(node.cuts_ref, CutScore(length(node.cutpool.hashMap[k]), viol, k))
                                con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                            end
                        end
                    end
                end
                if applied 
                    l = r+1 ; break
                end 
            end
        end
    end

end

"""
Compute and stock the relaxed bound set (i.e. the LP relaxation) of the (sub)problem defined by the given node.
Return `true` if the node is pruned by infeasibility.
"""
function LPRelaxByDicho(node::Node, pb::BO01Problem, round_results, verbose ; args...)
    start = time()

    assignment = getPartialAssign(node)
    setBounds(pb, assignment)


    pruned = compute_LBS(node, pb, round_results, verbose; args)
    pb.info.relaxation_time += (time() - start)

    if pruned 
        removeBounds(pb, assignment) ; return true 
    end

    if pb.param.cut_activated
        start_cuts = time()

        # add valid cuts constraints then re-optimize 
        start_processing = time()
        loadingCutInPool( node, pb)         # complexity O(pt ⋅ cuts)
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

        start_dicho = time()
        pruned = compute_LBS(node, pb, round_results, verbose; args)
        pb.info.cuts_infos.times_calling_dicho += (time() - start_dicho)

        if length(node.RBS.natural_order_vect) > 1
            pruned = MP_cutting_planes(node, pb, round_results, verbose ; args...)

        elseif length(node.RBS.natural_order_vect) == 1
            pb.info.cuts_infos.ite_total += 1 
            (new_x, cut_found) = SP_cut_off(1, node, pb, round_results, verbose ; args...)
            if cut_found && new_x != node.RBS.natural_order_vect.sols[1].xEquiv[1]
                node.RBS.natural_order_vect.sols[1].xEquiv[1] = new_x[:]
                node.RBS.natural_order_vect.sols[1].y = [pb.c[1, 1] + pb.c[1, 2:end]'*new_x , pb.c[2, 1] + pb.c[2, 2:end]'*new_x]
                node.RBS.natural_order_vect.sols[1].is_binary = isBinary(new_x)
            end
        end

        # retrieve applied valid cuts 
        start_processing = time()
        for con in node.con_cuts
            if JuMP.is_valid( pb.m, con)
                JuMP.delete( pb.m, con) ; JuMP.unregister( pb.m, :con) # remove the symbolic reference
            end
        end
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

        pb.info.cuts_infos.times_total_for_cuts += (time() - start_cuts)

        # println("-----------------------------------")
        # @info "node $(node.num) cutpool !"
        # for (k,v) in node.cutpool.hashMap
        #     println(k, " => ", size(v, 1))
        # end
        # @info " node $(node.num) cuts : $(length(node.con_cuts))"
        # println("-----------------------------------")
    end

    removeBounds(pb, assignment)

    return pruned
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