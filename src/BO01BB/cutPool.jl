# This file contains functions related to cutPool structure.

# mutable struct CutScore
#     cut_ind::Int64          # index of the referred cut in the cut pool 
#     viol::Float64           # the amount of the violation of the referred cut 
#     k::Int64
# end

# function CutScore()
#     return CutScore(0, 0.0, 0)
# end

# function push_cutScore(cut_refs::Vector{CutScore}, cut_score::CutScore)
#     if length(cut_refs) == 0
#         push!(cut_refs, cut_score) ; return true
#     end

#     for cs in cut_refs 
#         if cs.k==cut_score.k && cs.cut_ind == cut_score.cut_ind return false end 
#     end
#     push!(cut_refs, cut_score) ; return true 
# end


function hashing_ineq(x::Vector{Int64})
    # if x[1]==0 return 0 end 
    # return round(Int64, sum([(i-1)*x[i] for i=2:length(x)])/x[1])
    n = length(x)
    return (sum([(i-1)*x[i] for i=2:n]) + x[1])%n
end


mutable struct Cut
    row::Vector{Int64}
    hash_k::Int64
    sum_coeff::Int64
    len_coeff_nonul::Int64
    min_val_indx::Int64
    min_val::Int64
    max_val_indx::Int64
    max_val::Int64
    acc_even_coeff::Int64
    acc_odd_coeff::Int64
end

function Cut()
    return Cut(Vector{Int64}(), 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

function Cut(x::Vector{Int64})
    c = Cut() ; c.row = x
    c.hash_k = hashing_ineq(x) ; n=length(x)
    c.sum_coeff = 0 ; c.len_coeff_nonul = 0
    c.min_val = x[1] ; c.min_val_indx = 0
    c.max_val = 0 ; c.max_val_indx = 0
    c.acc_even_coeff = 0 ; c.acc_odd_coeff = 0

    for i=2:n 
        if x[i] != 0
            c.sum_coeff += x[i] ; c.len_coeff_nonul += 1

            if x[i] > c.max_val
                c.max_val = x[i] ; c.max_val_indx = i 
            end

            if x[i] < c.min_val
                c.min_val = x[i] ; c.min_val_indx = i 
            end

            if i%2==1
                c.acc_odd_coeff += (i-1)*x[i] 
            else
                c.acc_even_coeff += (i-1)*x[i] 
            end
        end
    end
    return c
end

function Base.isequal(c1::Cut, c2::Cut)
    @assert length(c1.row)> 0 && length(c2.row) > 0
    if c1.hash_k != c2.hash_k return false end 
    if c1.sum_coeff != c2.sum_coeff || c1.len_coeff_nonul != c2.len_coeff_nonul return false end 
    if c1.min_val != c2.min_val || c1.min_val_indx != c2.min_val_indx return false end 
    if c1.max_val != c2.max_val || c1.max_val_indx != c2.max_val_indx return false end 
    if c1.acc_even_coeff != c2.acc_even_coeff || c1.acc_odd_coeff != c2.acc_odd_coeff return false end 
    return true
end


mutable struct CutPool
    hashMap::Dict{Int64, Vector{Cut}}
end

function CutPool()
    return CutPool(Dict{Int64, Vector{Cut}}())
end

"""
Push the given `cut` in the cut pool only if `cut` is not redundent. 
Return `true` if `cut` is successfully added.
"""
function Base.push!(cpool::CutPool, cut::Cut)    
    @assert length(cut.row)> 0
    k = cut.hash_k
    if !haskey(cpool.hashMap, k)
        cpool.hashMap[k] = Vector{Cut}()
        push!(cpool.hashMap[k], cut) ; return true
    end

    for c in cpool.hashMap[k]
        if isequal(c, cut) return false end
    end

    push!(cpool.hashMap[k], cut) ; return true
end

function total_cuts(cpool::CutPool)
    count = 0
    for (k, v) in cpool.hashMap
        count += length(v)
    end
    return count
end