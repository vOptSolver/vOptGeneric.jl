# This file contains functions related to cutPool structure.

mutable struct CutScore
    cut_ind::Int64          # index of the referred cut in the cut pool 
    viol::Float64           # the amount of the violation of the referred cut 
    age::Int64              # the number of iterations applied
end

function CutScore()
    return CutScore(0, 0.0, 0)
end

function push_cutScore(cut_refs::Vector{CutScore}, cut_score::CutScore)
    if length(cut_refs) == 0
        push!(cut_refs, cut_score) ; return true
    end

    for cs in cut_refs 
        if cs.cut_ind == cut_score.cut_ind return false end 
    end
    push!(cut_refs, cut_score) ; return true 
end


mutable struct CutPool
    tab::Vector{Vector{Int64}}
end

function CutPool()
    return CutPool(Vector{Vector{Int64}}())
end

"""
Push the given `cut` in the cut pool only if `cut` is not redundent. 
Return `true` if `cut` is successfully added.
"""
function Base.push!(cpool::CutPool, cut::Vector{Int64})
    if size(cpool.tab, 1) == 0
        push!(cpool.tab, cut) ; return true
    end

    for c in cpool.tab
        notequal = false
        for i=1:length(c)
            if !isapprox(c[i], cut[i]; atol=10^(-4))
                notequal = true ; break
            end
        end

        if !notequal return false end
    end

    push!(cpool.tab, cut)
    # @info "cut : ", cut
    # println("number of cuts : ", length(cpool.tab))
    # for cut in cpool.tab
    #     println(cut)
    # end
    return true
end