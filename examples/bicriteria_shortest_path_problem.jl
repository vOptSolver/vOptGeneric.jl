using vOptGeneric
using GLPK,GLPKMathProgInterface
m = vModel(solver = GLPKSolverMIP())
# using CPLEX
# m = vModel(solver = CplexSolver())

@variable(m, x[1:6,1:6], Bin)

M = 50
p1 =   [M 4 5 M M M;
        M M 2 1 2 7;
        M M M 5 2 M;
        M M 5 M M 3;
        M M M M M 4;
        M M M M M M]
p2 =   [M 3 1 M M M;
        M M 1 4 2 2;
        M M M 1 7 M;
        M M 1 M M 2;
        M M M M M 2;
        M M M M M M]


@addobjective(m, Min , sum(x[i,j]*p1[i,j] for i=1:6, j=1:6))
@addobjective(m, Min , sum(x[i,j]*p2[i,j] for i=1:6, j=1:6))

@constraint(m, eps, sum(x[1,j] for j = 1:6) == 1)
@constraint(m, sum(x[i,6] for i = 1:6) == 1)
@constraint(m, cstr[i=2:5], sum(x[i,j] for j = 1:6) - sum(x[j,i] for j =1:6) == 0)

solve(m, method=:epsilon, step=0.5)

#get the results and print/plot them
Y_N = getY_N(m)

for n = 1:length(Y_N)
    X = getvalue(x, n)
    for ind in find(X)
        i,j = ind2sub(X,ind)
        print("$i->$j ")
    end
    println("| z = ",Y_N[n])
end

# using PyPlot
# f1 = map(y -> y[1], Y_N)
# f2 = map(y -> y[2], Y_N)
# xlabel("z1")
# ylabel("z2")
# plot(f1,f2,"bx", markersize = "6")
# !isinteractive() && show()
