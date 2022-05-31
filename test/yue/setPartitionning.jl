# Bi-objective set partitionning problem (2-SPA)


# ---- Packages to use
using JuMP, GLPK

include("../../src/vOptGeneric.jl")
using .vOptGeneric


# ---- Parser reading an instance of 2-SPA (format of instances compliant with vOptLib)
function loadInstance2SPA(fname::String)
    println("reading ", fname)
    f = open(fname)
    nbctr, nbvar = parse.(Int, split(readline(f))) # number of constraints , number of variables
    A = zeros(Int, nbctr, nbvar)                   # matrices of constraints
    c1 = zeros(Int, nbvar)                         # 1st vector of costs
    c2 = zeros(Int, nbvar)                         # 2nd vector of costs
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c1[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                c2[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                A[j,i] = 1
            end
        end
    end
    close(f)
    return c1, c2, A
end


# ---- Compute the set of non-dominated points Y_N of a 2SPA with vOptGeneric
function computeYNfor2SPA(  nbvar::Int,
                            nbctr::Int,
                            A::Array{Int,2},
                            c1::Array{Int,1},
                            c2::Array{Int,1}
                         )

    # ---- setting the model
    model = vModel( GLPK.Optimizer )
    JuMP.set_silent(model)

    @variable(model, x[1:nbvar], Bin)
    @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    @addobjective(model, Min, sum(c1[i]*x[i] for i in 1:nbvar))
    @addobjective(model, Min, sum(c2[i]*x[i] for i in 1:nbvar))

    # ---- Invoking the solver (branch and bound method)
    vSolve( model, method=:bb, verbose=true )

    # # ---- Invoking the solver (epsilon constraint method)
    # vSolve(model, method=:epsilon, step = 0.5, verbose=true)

    # # ---- Querying the results
    # Y_N = getY_N(model)
    # for i = 1:length(Y_N)
    #     X = value.(x, i)
    #     print("X = ", findall(elt -> elt ≈ 1, X))
    #     println(" | Z = ",Y_N[i])
    # end
end


# ---- Entry point
function main(fname::String)

    print("Compute Y_N with vOPtGeneric for a set partitionning problem with 2 objectives (2-SPA) \n\n")

    # load a numerical instance of 2SPA ----------------------------------------
    c1, c2, A = loadInstance2SPA(fname)
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    # compute the set of non-dominated points Y_N ------------------------------
    computeYNfor2SPA(nbvar, nbctr, A, c1, c2)
end

# ---- Example
main("biodidactic5.txt")

# method=:epsilon, step = 0.5,

# X = [3, 6, 7] | Z = [14.0, 48.0]
# X = [3, 4, 6] | Z = [16.0, 39.0]
# X = [1, 5, 7] | Z = [18.0, 16.0]
# X = [1, 4, 5] | Z = [20.0, 7.0]
# X = [5, 8] | Z = [21.0, 4.0]

# method=:bb

# incumbent : NaturalOrderVector[
# Solution( 
#  xEquiv = [[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]]
#  y = [21.0, 4.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [20.0, 7.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0]]
#  y = [18.0, 16.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0]]
#  y = [14.0, 48.0]
#  is_binary ? true )
# ]

