using vOptGeneric
using Base.Test
using GLPK, GLPKMathProgInterface


# write your own tests here
m = MultiModel(solver = GLPKSolverMIP())

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

@constraint(m, sum(x[1,j] for j = 1:6) == 1)
@constraint(m, sum(x[i,6] for i = 1:6) == 1)
@constraint(m, cstr[i=2:5], sum(x[i,j] for j = 1:6) - sum(x[j,i] for j =1:6) == 0)

solve(m, method=:eps, ϵ=0.1)

#get the results and print/plot them
md = getMultiData(m)
Y_N = md.Y_N
X_E = md.X_E

f1 = map(x -> x[1], Y_N)
f2 = map(x -> x[2], Y_N)

@test f1 == [8.0,10.0,11.0,13.0]
@test f2 == [9.0,7.0,5.0,4.0]