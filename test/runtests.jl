using vOptGeneric
using Test, LinearAlgebra
using GLPK, GLPKMathProgInterface

function run_tests()
        m = vModel(solver = GLPKSolverMIP())

        @variable(m, x[1:6,1:6], Bin)
        @variable(m, y, Bin)
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

        solve(m, method=:epsilon)

        Y_N = getY_N(m)
        f1 = map(x -> x[1], Y_N)
        f2 = map(x -> x[2], Y_N)

        @test f1 == [8.0,10.0,11.0,13.0]
        @test f2 == [9.0,7.0,5.0,4.0]
        @test getvalue(x, 1) == [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
        printX_E(m)

        writeMOP(m, "test.txt")
        m = parseMOP("test.txt", solver = GLPKSolverMIP())
        solve(m, method=:dichotomy)
        Y_N = getY_N(m)
        f1 = map(x -> x[1], Y_N)
        f2 = map(x -> x[2], Y_N)
        @test f1 == [8.0, 11.0, 13.0]
        @test f2 == [9.0, 5.0, 4.0]

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

        solve(m, method=:chalmet)

        Y_N = getY_N(m)
        f1 = map(x -> x[1], Y_N)
        f2 = map(x -> x[2], Y_N)

        @test f1 == [918.0, 924.0, 927.0, 934.0, 935.0, 940.0, 943.0, 948.0, 949.0, 952.0, 955.0]
        @test f2 == [984.0, 975.0, 972.0, 971.0, 947.0, 943.0, 940.0, 939.0, 915.0, 909.0, 906.0]
        @test getvalue(x, 1) == [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0]


        solve(m, method=:lex)
        Y_N = getY_N(m)
        @test Y_N == [[955.0, 906.0], [918.0, 984.0]]
end

run_tests()