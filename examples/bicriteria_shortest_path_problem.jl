using vOptGeneric, JuMP
using Cbc

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

m = vModel(with_optimizer(Cbc.Optimizer, logLevel=0))
@variable(m, x[1:6,1:6], Bin)

@addobjective(m, Min , sum(x[i,j]*p1[i,j] for i=1:6, j=1:6))
@addobjective(m, Min , sum(x[i,j]*p2[i,j] for i=1:6, j=1:6))
@constraint(m, eps, sum(x[1,j] for j = 1:6) == 1)
@constraint(m, sum(x[i,6] for i = 1:6) == 1)
@constraint(m, cstr[i=2:5], sum(x[i,j] for j = 1:6) - sum(x[j,i] for j =1:6) == 0)

vSolve(m, method=:epsilon, step=0.5, verbose=true)

#get the results and print them
Y_N = getY_N(m)
for n = 1:length(Y_N)
    X = value.(x, n)
    for ind in findall(val -> val ≈ 1, X)
        i,j = ind.I
        print("$i->$j ")
    end
    println("| z = ",Y_N[n])
end