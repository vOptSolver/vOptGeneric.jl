using vOptGeneric
using GLPK,GLPKMathProgInterface
m = vModel(solver = GLPKSolverMIP())
# using CPLEX
# m = vModel(solver = CplexSolver())

p1 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62]
p2 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
w = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99]
c = 900
size = length(p1)

@variable(m, x[1:size], Bin)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2))
@constraint(m, dot(x, w) <= c)

solve(m, method=:dichotomy)

Y_N = getY_N(m)
for n = 1:length(Y_N)
    X = getvalue(x, n)
    print(find(X))
    println("| z = ",Y_N[n])
end

# using PyPlot
# f1 = map(y -> y[1], Y_N)
# f2 = map(y -> y[2], Y_N)
# xlabel("z1")
# ylabel("z2")
# plot(f1,f2,"bx", markersize = "6")
# !isinteractive() && show()