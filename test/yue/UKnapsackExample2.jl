# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Example 2 with 20 variables


# ---- Packages to use
using JuMP, GLPK, LinearAlgebra

include("../../src/vOptGeneric.jl")
using .vOptGeneric

# ---- Values of the instance to solve
p1 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62]
p2 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
w  = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99]
c  = 900
size = length(p1)


# ---- setting the model
m = vModel( GLPK.Optimizer ) ; JuMP.set_silent(m)
@variable( m, x[1:size], Bin )
@addobjective( m, Max, dot(x, p1) )
@addobjective( m, Max, dot(x, p2) )
@constraint( m, dot(x, w) <= c )

# ---- Invoking the solver (branch and bound method)
infos = vSolve( m, method=:bb, verbose=true )
println("infos : ", infos)

# # ---- Invoking the solver (epsilon constraint method)
# vSolve( m, method=:epsilon, step=0.5, verbose=true )


# ---- Querying the results
Y_N = getY_N(m)
println("number Y_N = ", length(Y_N))


# # ---- Displaying the results (X_E and Y_N)
# for n = 1:length(Y_N)
#     X = value.(x, n)
#     print(findall(elt -> elt ≈ 1, X))
#     println("| z = ",Y_N[n])
# end

# method=:epsilon, step=0.5; Output : 

# [2, 3, 4, 5, 6, 8, 9, 12, 15, 16, 18, 19]| z = [918.0, 984.0]
# [2, 3, 5, 6, 8, 10, 11, 12, 16, 17, 19]| z = [924.0, 975.0]
# [2, 3, 5, 6, 8, 9, 10, 11, 12, 16, 17]| z = [927.0, 972.0]
# [2, 3, 5, 6, 8, 10, 11, 12, 15, 16, 17]| z = [934.0, 971.0]
# [2, 5, 6, 8, 9, 10, 11, 12, 15, 16, 17]| z = [935.0, 947.0]
# [2, 3, 5, 6, 8, 10, 11, 15, 16, 17, 19]| z = [940.0, 943.0]
# [2, 3, 5, 6, 8, 9, 10, 11, 15, 16, 17]| z = [943.0, 940.0]
# [1, 2, 3, 5, 6, 8, 10, 11, 15, 16, 17]| z = [948.0, 939.0]
# [1, 2, 5, 6, 8, 9, 10, 11, 15, 16, 17]| z = [949.0, 915.0]
# [2, 3, 5, 6, 10, 11, 14, 15, 16, 17, 19]| z = [952.0, 909.0]
# [2, 3, 5, 6, 9, 10, 11, 14, 15, 16, 17]| z = [955.0, 906.0]

# method=:bb ; Output : 
# incumbent : NaturalOrderVector[
# Solution( 
#  xEquiv = [[0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0]]
#  y = [924.0, 975.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0]]
#  y = [918.0, 984.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [927.0, 972.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [934.0, 971.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [920.0, 936.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [934.0, 935.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [933.0, 925.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [934.0, 920.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]]
#  y = [934.0, 904.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]]
#  y = [936.0, 900.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0]]
#  y = [942.0, 886.0]
#  is_binary ? true )
# ]



# using PyPlot
# f1 = map(y -> y[1], Y_N)
# f2 = map(y -> y[2], Y_N)
# xlabel("z1")
# ylabel("z2")
# plot(f1,f2,"bx", markersize = "6")
# !isinteractive() && show()
