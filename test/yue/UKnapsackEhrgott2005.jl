# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Exercise 10.2 page 290 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using JuMP, GLPK

include("../../src/vOptGeneric.jl")
using .vOptGeneric

# ---- Values of the instance to solve
p1 = [10, 3,  6, 8, 2]  # coefficients's vector of the objective 1
p2 = [12, 9, 11, 5, 6]  # coefficients's vector of the objective 2
w  = [ 4, 5,  2, 5, 6]  # coefficients's vector of weights
c  = 17                 # nominal capacity
n  = length( p1 )       # number of items


# ---- setting the model
biukp = vModel( GLPK.Optimizer ) #; JuMP.set_silent( biukp )
@variable( biukp, x[1:n], Bin )
@addobjective( biukp, Max, sum( p1[j]*x[j] for j=1:n ) )
@addobjective( biukp, Max, sum( p2[j]*x[j] for j=1:n ) )
@constraint( biukp, sum( w[j]*x[j] for j=1:n ) <= c )


# ---- Invoking the solver (branch and bound method)
vSolve( biukp, method=:bb, verbose=true )

# # ---- Invoking the solver (dichotomy method)
# vSolve( biukp, method=:dicho, verbose=true )

# # ---- Querying the results
# Y_N = getY_N( biukp )


# # ---- Displaying the results (X_{SE} and Y_{SN})
# for i = 1:length(Y_N)
#     X = value.(x, i)
#     print("X = ", findall(elt -> elt ≈ 1, X))
#     println(" | Z = ",Y_N[i])
# end



# method=:dicho, verbose=true

# Output : 

# solving for z1
# solving for z2
# solving for 1.0*f1 + 6.0*f2
# X = [1, 2, 3, 5] | Z = [21.0, 38.0]
# X = [1, 2, 3, 4] | Z = [27.0, 37.0]


# method=:bb, verbose=true

# Output : 

# incumbent : NaturalOrderVector[
# Solution( 
#  xEquiv = [[1.0, 1.0, 1.0, 0.0, 1.0]]
#  y = [21.0, 38.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 1.0, 1.0, 1.0, 0.0]]
#  y = [27.0, 37.0]
#  is_binary ? true )
# ]

