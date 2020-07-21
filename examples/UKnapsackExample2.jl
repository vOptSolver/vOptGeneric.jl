# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Example 2 with 20 variables


# ---- Packages to use
using vOptGeneric, JuMP, GLPK, LinearAlgebra


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


# ---- Invoking the solver (epsilon constraint method)
vSolve( m, method=:epsilon, step=0.5, verbose=true )


# ---- Querying the results
Y_N = getY_N(m)


# ---- Displaying the results (X_E and Y_N)
for n = 1:length(Y_N)
    X = value.(x, n)
    print(findall(elt -> elt ≈ 1, X))
    println("| z = ",Y_N[n])
end


# using PyPlot
# f1 = map(y -> y[1], Y_N)
# f2 = map(y -> y[2], Y_N)
# xlabel("z1")
# ylabel("z2")
# plot(f1,f2,"bx", markersize = "6")
# !isinteractive() && show()
