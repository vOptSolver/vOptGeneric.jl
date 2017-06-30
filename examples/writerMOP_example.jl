using vOptGeneric
m = vModel()

@variable(m, x[1:6,1:6], Bin)
@variable(m, 0 <= y <= 10, Int)

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

@addobjective(m, Min , y + sum(x[i,j]*p1[i,j] for i=1:6, j=1:6))
@addobjective(m, Min , y + sum(x[i,j]*p2[i,j] for i=1:6, j=1:6))

@constraint(m, sum(x[1,j] for j = 1:6) == 1)
@constraint(m, sum(x[i,6] for i = 1:6) == 1)
@constraint(m, cstr[i=2:5], sum(x[i,j] for j = 1:6) - sum(x[j,i] for j =1:6) == 0)
@constraint(m, y <= 5)

writeMOP(m, "bicriteria_shortest_path.mop")