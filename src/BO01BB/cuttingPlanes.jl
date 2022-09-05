# This file contains functions of cutting planes algorithm.
include("BBtree.jl")
include("../algorithms.jl")
include("struct.jl")
include("separators.jl")
include("cutPool.jl")

using JuMP, CPLEX 

const max_step = 2
const loop_limit = 5


"""
Compute the lower bound set of the LP polyhedron by dichotomy method.

Return `true` if this node is fathomed by infeasibility.
"""
function compute_LBS(node::Node, pb::BO01Problem, round_results, verbose ; args...)
    #------------------------------------------------------------------------------
    # solve the LP relaxation by dichotomy method including the partial assignment
    #------------------------------------------------------------------------------
    solve_dicho(pb.m, round_results, false ; args...)
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
        # pb.info.status = MOI.INFEASIBLE
        return true
    end

    # construct/complete the relaxed bound set
    node.RBS = RelaxedBoundSet() #; node.objs = Vector{JuMP.GenericAffExpr}()
    for i = 1:length(vd_LP.Y_N)
        push!(node.RBS.natural_order_vect, Solution(vd_LP.X_E[i], vd_LP.Y_N[i]))
    end
    for i=1:length(node.RBS.natural_order_vect)-1
        node.RBS.segments[node.RBS.natural_order_vect.sols[i]] = true
    end

    return false
end


"""
A single point cut off algorithm. 

Returns :
    - the fractional/integer point
    - true, if new cut(s) has benn found 
"""
function SP_cut_off(i::Int64, node::Node, pb::BO01Problem, round_results, verbose ; args...)
    x_star = node.RBS.natural_order_vect.sols[i].xEquiv[1]
    if node.RBS.natural_order_vect.sols[i].is_binary return (x_star, false) end 

    start_sep = time()
    cuts = SP_KP_heurSeparator2(x_star, pb.A, pb.b)
    pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

    if length(cuts) > 0
        for cut in cuts
            start_pool = time()
            ineq = Cut(cut)
            if push!(node.cutpool, ineq)
                pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
                con = JuMP.@constraint(pb.m, cut[2:end]'*pb.varArray ≤ cut[1]) ; push!(node.con_cuts, con)
            end
            pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

        end
        return (x_star, true)
    end

    # # call generator
    # start_sep = time()
    # (isValidCut, α, _) = SP_CG_separator(x_star, pb.A, pb.b)
    # pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

    # if isValidCut
    #     # @info " ------------------------- cut found"
    #     start_pool = time()
    #     ineq = Cut(α)
    #     if push!(node.cutpool, ineq)
    #         pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
    #         con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
    #     end
    #     pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)
    #     return (x_star, true)
    # end

    return (x_star, false)
end


"""
Cutting planes scheme for the multi-point cuts. For now we assume that every point `y` in criteria space has only
one corresponding vector `x` in decision space. 

Return ture if the node is infeasible after adding cuts.
"""
function MP_cutting_planes(node::Node, pb::BO01Problem, round_results, verbose ; args...)
    numVars = length(pb.varArray) ; numRows = size(pb.A, 1)
    LBS = node.RBS.natural_order_vect.sols

    ite = 0
    while ite < loop_limit 
        ite += 1 ; pb.info.cuts_infos.ite_total += 1 
        
        # @info "MP_cutting_planes ite = $ite"
        # ------------------------------------------------------------------------------
        # 1. generate multi-point cuts if has any, or single-point cut off
        # ------------------------------------------------------------------------------
        l = 1 ; cut_counter = 0

        while l ≤ length(LBS)
            if LBS[l].is_binary 
                l += 1 ; continue    
            end

            for ∇ = max_step:-1:0 
                if ∇ == 0
                    # @info "MP_cutting_planes calling `SP_cut_off`"
                    (_, new_cut) = SP_cut_off(l, node, pb, round_results, verbose ; args...) 
                    if new_cut cut_counter += 1 end 
                    l += 1
                else 
                    r = l+∇
                    if r > length(LBS) || LBS[r].is_binary continue end

                    start_sep = time()
                    cuts = MP_KP_heurSeparator2(LBS[l].xEquiv[1], LBS[r].xEquiv[1], pb.A, pb.b)
                    pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

                    if length(cuts) > 0
                        cut_counter += (∇+1)
                        for cut in cuts
                            start_pool = time()
                            ineq = Cut(cut)
                            if push!(node.cutpool, ineq)
                                pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                                con = JuMP.@constraint(pb.m, cut[2:end]'*pb.varArray ≤ cut[1]) ; push!(node.con_cuts, con)
                            end
                            pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

                        end
                        l = r + 1 ; break
                    end

                    # start_sep = time()
                    # (isValidCut, α, _) = MP_CG_separator(LBS[l].xEquiv[1], LBS[r].xEquiv[1], pb.A, pb.b)
                    # pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

                    # if isValidCut
                    #     cut_counter += (∇+1)
                    #     # @info " ---------------------------- cut found "
                    #     start_pool = time()
                    #     ineq = Cut(α)
                    #     if push!(node.cutpool, ineq)
                    #         pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                    #         con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                    #     end
                    #     pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)
                    #     l = r + 1 ; break
                    # end

                end
    
            end
        end

        # --------------------------------------------------
        # 2. stop if no more valid cut can be found
        # --------------------------------------------------
        if cut_counter/length(LBS) < 0.4
            return false 
        end

        # ---------------------------------------------------
        # 3. otherwise, re-optimize by solving dicho -> LBS
        # ---------------------------------------------------
        start_dicho = time()
        pruned = compute_LBS(node, pb, round_results, verbose; args)
        pb.info.cuts_infos.times_calling_dicho += (time() - start_dicho)
        LBS = node.RBS.natural_order_vect.sols

        # in case of infeasibility
        if pruned return true end

        # in case of integrity
        all_integers = true
        for sol in node.RBS.natural_order_vect.sols
            if !sol.is_binary 
                all_integers = false ; break 
            end
        end
        if all_integers return false end

    end # loop while

    return false 
end
