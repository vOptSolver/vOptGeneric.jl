# This file contains functions related to node fathoming.

include("cuttingPlanes.jl")


"""
Given a point `x_star`, iterate all valid cuts of parent node and stock the references 
in the current node.
"""
function loadingCutInPool(node::Node, pb::BO01Problem)
    if isRoot(node) return end 

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
                    for cut in cuts
                        α = cut.row
                        violationₗ = maximum([ (xₗ_star'*α[2:end] - α[1]), 0.0 ])
                        if violationₗ > 0.0
                            # ineq = Cut(α)
                            if push!(node.cutpool, cut)# && push_cutScore(node.cuts_ref, CutScore(length(node.cutpool.hashMap[k]), violationₗ, k))
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
                    for cut in cuts 
                        α = cut.row
                        violationₗ = maximum([ (xₗ_star'*α[2:end] - α[1]), 0.0 ])
                        violationᵣ = maximum([ (xᵣ_star'*α[2:end] - α[1]), 0.0 ])
                        viol = maximum([violationₗ, violationᵣ])
                        if viol > 0.0
                            applied = true
                            # ineq = Cut(α)
                            if push!(node.cutpool, cut)# && push_cutScore(node.cuts_ref, CutScore(length(node.cutpool.hashMap[k]), viol, k))
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
    objcons = setVarObjBounds(node, pb) ; num_var = length(pb.varArray)

    # ------------------------
    # apply valid cuts 
    # ------------------------
    if pb.param.cut_activated #&& node.depth < num_var/3
        start_cuts = time() ; pruned = false 

        # step 1 : calculate LBS of the actual sub-problem
        start = time()
        pruned = compute_LBS(node, pb, round_results, verbose; args)
        pb.info.relaxation_time += (time() - start)

        if pruned 
            removeVarObjBounds(node, pb, objcons) ; return true 
        end

        @assert length(node.RBS.natural_order_vect) > 0 "valid LBS is empty !"

        # step 2 : add valid cuts constraints then re-optimize 
        start_processing = time()
        loadingCutInPool( node, pb)         # complexity O(pt ⋅ cuts)
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

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

        # step 3 : retrieve applied valid cuts 
        start_processing = time()
        for con in node.con_cuts
            if JuMP.is_valid( pb.m, con)
                JuMP.delete( pb.m, con) ; JuMP.unregister( pb.m, :con) # remove the symbolic reference
            end
        end
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

        pb.info.cuts_infos.times_total_for_cuts += (time() - start_cuts)

        removeVarObjBounds(node, pb, objcons) ; return pruned
    else
        start = time()
        pruned = compute_LBS(node, pb, round_results, verbose; args)
        pb.info.relaxation_time += (time() - start)

        removeVarObjBounds(node, pb, objcons) ; return pruned
    end

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
        pb.info.update_incumb_time += (time() - start) ; return true
    end
    pb.info.update_incumb_time += (time() - start) ; return false
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
function fullyExplicitDominanceTest(node::Node, global_incumbent::IncumbentSet)
    @assert length(node.RBS.natural_order_vect) > 0 "relaxed bound set is empty for node $(node.num)"

    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(global_incumbent.natural_order_vect) == 0 return false end

    # if node.EPB     # consider a "local" upper bound sets 
    #     incumbent = IncumbentSet()
    #     for u in global_incumbent.natural_order_vect.sols
    #         if u.y[1] ≤ node.nadirPt[1] && u.y[2] ≤ node.nadirPt[2]
    #             push!(incumbent.natural_order_vect, u)
    #         end
    #     end
    # else 
    #     incumbent = global_incumbent
    # end

    incumbent = global_incumbent

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
    # two extreme points of LBS
    ptl = node.RBS.natural_order_vect.sols[1] ; ptr = node.RBS.natural_order_vect.sols[end]

    # Case 1 :  if only one feasible point in UBS 
    if length(incumbent.natural_order_vect) == 1 
        # who dominates the ideal point 
        if incumbent.natural_order_vect.sols[1].y[1] ≤ ptr.y[1] && incumbent.natural_order_vect.sols[1].y[2] ≤ ptl.y[2]
            return true
        else
            # Pareto branching 
            #TODO : check => Pareto branching ... 
            return false 
        end
        # return false 
    end

    # Case 2 : otherwise, do the pairwise comparison of the local nadir points with LBS  
    (nadir_pts, fathomed) = getNadirPoints(incumbent, ptl, ptr)
    if fathomed return true end

    # test range condition necessary 1 : LBS ⊆ UBS 
    u_l = incumbent.natural_order_vect.sols[1] ; u_r = incumbent.natural_order_vect.sols[end]

    sufficient = (u_l.y[2] < ptl.y[2] && u_r.y[1] < ptr.y[1])

    if !sufficient return false end

    # test condition necessary 2 : LBS ≤/dominates UBS 
    fathomed = true
    # iterate of all local nadir points
    for u ∈ nadir_pts.sols
        existence = false ; compared = false

        # case 1 : if u is dominates the ideal point of LBS 
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        end

        # # case 2 : if u is worse than the "worst nadir point" of LBS # todo: actually doesn't work  
        # if u.y[1] > ptl.y[1] && u.y[2] > ptr.y[2]
        #     return false
        # end

        # case 3 : complete pairwise comparison
        for i=1:length(node.RBS.natural_order_vect)-1              # ∀ segment l ∈ LBS 

            sol_l = node.RBS.natural_order_vect.sols[i] ; sol_r = node.RBS.natural_order_vect.sols[i+1]

            if (u.y[1] > sol_l.y[1] || u.y[1] < sol_r.y[1]) && (u.y[2] > sol_r.y[2] || u.y[2] < sol_l.y[2])
                continue
            end
            
            λ = [sol_r.y[2] - sol_l.y[2], sol_l.y[1] - sol_r.y[1]]      # normal to the segment

            compared = true

            if λ'*u.y < λ'*sol_r.y #&& λ'*u.y < λ'*sol_l.y
                existence = true ; break
            end
        end
        
        # case 4 : condition dominance violated, then stock the non-dominated local nadir pts to prepare EPB
        if compared && !existence 
            fathomed = false

            if !isRoot(node) && (u.y in node.pred.localNadirPts || u.y == node.pred.nadirPt || u.y == node.nadirPt)    # the current local nadir pt is already branched 
            # if u.y in node.pred.localNadirPts
                node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed
            else#if EPB_decider(u.y, ptl, ptr)
                push!(node.localNadirPts, u.y)
            # else
            #     node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed
            end
        end

        if !compared && (u.y[1] ≥ ptr.y[1] && u.y[2] ≥ ptl.y[2])
            node.localNadirPts = Vector{Vector{Float64}}()              # no need to (extended) pareto branching
            return false
        end
    end

    return fathomed
end


"""
Return the orthogonal distance between a point `p` and a segment defined by two points `extl` and `extr`.
"""
function orthogonal_dist(p::Vector{Float64}, extl::Vector{Float64}, extr::Vector{Float64})
    return abs((extr[1] - extl[1]) * (extl[2] - p[2]) - (extl[1] - p[1]) * (extr[2] - extl[2])) / sqrt((extr[1] - extl[1])^2 + (extr[2] - extl[2])^2) 
end


"""
Given a non-dominated nadir point, return `true` if the decider EP - branch on it, considering the average distance of nadir point from LBS.
"""
function EPB_decider( node::Node)
    ptl = node.RBS.natural_order_vect.sols[1].y ; ptr = node.RBS.natural_order_vect.sols[end].y
    worst_nadir_pt = [ptl[1], ptr[2]] ; dist_LBS = orthogonal_dist(worst_nadir_pt, ptl, ptr)
    
    avg_dist = 0.0
    for nadir_pt in node.localNadirPts
        dist_nadir = orthogonal_dist(nadir_pt, ptl, ptr)
        avg_dist += dist_nadir/dist_LBS
    end
    
    return (avg_dist/length(node.localNadirPts)) ≤ 1/2
end