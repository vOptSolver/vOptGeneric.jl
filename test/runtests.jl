using JuMP, vOptGeneric, Cbc, Combinatorics
using Test, LinearAlgebra

m = vModel(with_optimizer(Cbc.Optimizer, logLevel=0))
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

@addobjective(m, Min, sum(x[i,j]*p1[i,j] for i=1:6, j=1:6) + 5)
@addobjective(m, Min, sum(x[i,j]*p2[i,j] for i=1:6, j=1:6) - 5)
@constraint(m, sum(x[1,j] for j = 1:6) == 1)
@constraint(m, sum(x[i,6] for i = 1:6) == 1)
@constraint(m, cstr[i=2:5], sum(x[i,j] for j = 1:6) - sum(x[j,i] for j =1:6) == 0)

vSolve(m, method=:epsilon)

Y_N = getY_N(m)
f1, f2 = first.(Y_N), last.(Y_N)
@test f1 == [13, 15, 16, 18] && f2 == [4, 2, 0, -1]
@test value.(x, 1) == [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
printX_E(m)

######################

# to test once https://github.com/JuliaOpt/JuMP.jl/issues/1956 is fixed
@test_broken begin
    m_copy = copy(m)        
    # vSolve(m_copy, with_optimizer(Cbc.Optimizer, logLevel=0), method=:epsilon)
    # Y_N = getY_N(m)
    # f1 = map(x -> x[1], Y_N) ; f2 = map(x -> x[2], Y_N)
    # @test f1 == [8.0,10.0,11.0,13.0] && f2 == [9.0,7.0,5.0,4.0]
    # @test value.(x, 1) == [0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    # printX_E(m_copy)
end

######################

# tests DenseAxisArrays indexing
# http://www.juliaopt.org/JuMP.jl/v0.19.1/variables/#variable_jump_arrays-1

m = vModel(with_optimizer(Cbc.Optimizer, logLevel=0))
cities = ["Paris", "New-York", "Madrid"]
@variable(m, x[cities], Bin)
@constraint(m, constr, x["Paris"] >= 1)
@addobjective(m, Min, sum(x))
@addobjective(m, Max, sum(x))
vSolve(m, method=:epsilon)
@show value.(x, 1)
@test value.(x, 1)["Paris"] == 1.0
@test value.(x, 1)["New-York"] == 0.0
@test value.(x, 1)["Madrid"] == 0.0

######################

# TODO : write a test for SparseAxisArray
# http://www.juliaopt.org/JuMP.jl/v0.19.1/variables/#variable_sparseaxisarrays-1

######################

# TODO : re-implement wirteMOP and parseMOP
# writeMOP(m, "test.txt")
# m = parseMOP("test.txt", solver = GLPKSolverMIP())
# solve(m, method=:dichotomy)
# Y_N = getY_N(m)
# f1 = map(x -> x[1], Y_N)
# f2 = map(x -> x[2], Y_N)
# @test f1 == [8.0, 11.0, 13.0]
# @test f2 == [9.0, 5.0, 4.0]

######################

m = vModel(with_optimizer(Cbc.Optimizer, logLevel=0))

p1 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62]
p2 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
w = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99]
c = 900
size = length(p1)

@variable(m, x[1:size], Bin)
@addobjective(m, Max, dot(x, p1))
@addobjective(m, Max, dot(x, p2))
@constraint(m, dot(x, w) <= c)

vSolve(m, method=:dicho)

Y_N = getY_N(m)
f1, f2 = first.(Y_N), last.(Y_N)

@test f1 == [918.0, 934.0, 948.0, 955.0]
@test f2 == [984.0, 971.0, 939.0, 906.0]
@test value.(x, 1) ≈ [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0]

vSolve(m, method=:chalmet)

Y_N = getY_N(m)
f1, f2 = first.(Y_N), last.(Y_N)

@test f1 == [918., 924., 927., 934., 935., 940., 943., 948., 949., 952., 955.]
@test f2 == [984., 975., 972., 971., 947., 943., 940., 939., 915., 909., 906.]
@test value.(x, 1) ≈ [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0]

########################

m = vModel(with_optimizer(Cbc.Optimizer, logLevel=0))
@variable(m, 0 <= x[1:4] <= 2, Int)
@addobjective(m, Max, sum(x))
@addobjective(m, Min, sum(x))
@addobjective(m, Max, x[3])
@addobjective(m, Max, x[4])
@constraint(m, sum(x) <= 3)
vSolve(m, method=:lex)
Y_N = getY_N(m)

for (ind, perm) in enumerate(permutations(1:4, 4))
    if perm[1] == 2
        @test Y_N[ind] == [0., 0., 0., 0.]
    elseif perm[2] == 2 && perm[1] != 1
        @test maximum(Y_N[ind]) == 2.
    end
end