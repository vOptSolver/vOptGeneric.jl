# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Exercise 10.2 page 290 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using vOptGeneric, JuMP, GLPK


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


# ---- Invoking the solver (dichotomic method)
vSolve( biukp, method=:dichotomy, verbose=false )


# ---- Querying the results
Y_N = getY_N( biukp )


# ---- Displaying the results (X_{SE} and Y_{SN})
for i = 1:length(Y_N)
    X = value.(x, i)
    print("X = ", findall(elt -> elt ≈ 1, X))
    println(" | Z = ",Y_N[i])
end
