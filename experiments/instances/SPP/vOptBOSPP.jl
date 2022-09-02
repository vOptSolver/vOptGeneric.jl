
const verbose = false
const graphic = true

using JuMP, CPLEX
include("../../../src/vOptGeneric.jl")
using .vOptGeneric 

include("parserBOSPP.jl")


function solveDichotomy(fname::String)
    if !isfile(fname)
        @error "This file doesn't exist ! $fname"
    end

    result_dir = "../../results/SPP/" * split(fname, "/")[2]
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    folder = result_dir * "/" * string(:dicho)
    if !isdir(folder)
        mkdir(folder)
    end
    outputName = folder * "/" * split(fname, "/")[end]
    inst = readingBOSPP(fname)

    # ---- setting the model
    println("Building...")
    bospp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(bospp)
    @variable( bospp, x[1:inst.n], Bin )
    @addobjective( bospp, Max, sum(inst.C1[j] * x[j] for j=1:inst.n) )
    @addobjective( bospp, Max, sum(inst.C2[j] * x[j] for j=1:inst.n) )

    @constraint( bospp, cte[i=1:inst.m], sum(x[j] for j in inst.Cover[i]) <= 1)

    println("Solving...")

    start = time()
    vSolve( bospp, method=:dicho, verbose=false)
    total_time = round(time() - start, digits = 2)

    # ---- Querying the results
    println("Querying...")
    Y_N = getY_N( bospp )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( bospp )
    println("length X_E = ", length(X_E))

    # ---- Writing the results
    writeResults(inst.n, inst.m, fname, outputName, :dicho, Y_N, X_E; total_time)
end



function solve_epsilon(fname::String)
    if !isfile(fname)
        @error "This file doesn't exist ! $fname"
    end

    result_dir = "../../results/SPP/" * split(fname, "/")[2]
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    folder = result_dir * "/" * string(:epsilon)
    if !isdir(folder)
        mkdir(folder)
    end
    outputName = folder * "/" * split(fname, "/")[end]
    inst = readingBOSPP(fname)

    # ---- setting the model
    println("Building...")
    bospp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(bospp)
    @variable( bospp, x[1:inst.n], Bin )
    @addobjective( bospp, Max, sum(inst.C1[j] * x[j] for j=1:inst.n) )
    @addobjective( bospp, Max, sum(inst.C2[j] * x[j] for j=1:inst.n) )

    @constraint( bospp, cte[i=1:inst.m], sum(x[j] for j in inst.Cover[i]) <= 1)

    println("Solving...")

    start = time()
    vSolve( bospp, method=:epsilon, step=0.5, verbose=false)
    total_time = round(time() - start, digits = 2)

    # ---- Querying the results
    println("Querying...")
    Y_N = getY_N( bospp )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( bospp )
    println("length X_E = ", length(X_E))

    # ---- Writing the results
    writeResults(inst.n, inst.m, fname, outputName, :epsilon, Y_N, X_E; total_time)
end


function solveBOBB(fname::String)
    if !isfile(fname)
        @error "This file doesn't exist ! $fname"
    end

    result_dir = "../../results/SPP/" * split(fname, "/")[2]
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    folder = result_dir * "/" * string(:bb)
    if !isdir(folder)
        mkdir(folder)
    end
    outputName = folder * "/" * split(fname, "/")[end]
    inst = readingBOSPP(fname)

    # ---- setting the model
    println("Building...")
    bospp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(bospp)
    @variable( bospp, x[1:inst.n], Bin )
    @addobjective( bospp, Max, sum(inst.C1[j] * x[j] for j=1:inst.n) )
    @addobjective( bospp, Max, sum(inst.C2[j] * x[j] for j=1:inst.n) )

    @constraint( bospp, cte[i=1:inst.m], sum(x[j] for j in inst.Cover[i]) <= 1)

    println("Solving...")

    infos = vSolve( bospp, method=:bb, verbose=false )
    println(infos)

    # ---- Querying the results
    println("Querying...")
    Y_N = getY_N( bospp )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( bospp )
    println("length X_E = ", length(X_E))

    # ---- Writing the results
    writeResults(inst.n, inst.m, fname, outputName, :bb, Y_N, X_E; infos)
end



function solveBOBC(fname::String)
    if !isfile(fname)
        @error "This file doesn't exist ! $fname"
    end

    result_dir = "../../results/SPP/" * split(fname, "/")[2]
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    folder = result_dir * "/" * string(:bc)
    if !isdir(folder)
        mkdir(folder)
    end
    outputName = folder * "/" * split(fname, "/")[end]
    inst = readingBOSPP(fname)

    # ---- setting the model
    println("Building...")
    bospp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(bospp)
    @variable( bospp, x[1:inst.n], Bin )
    @addobjective( bospp, Max, sum(inst.C1[j] * x[j] for j=1:inst.n) )
    @addobjective( bospp, Max, sum(inst.C2[j] * x[j] for j=1:inst.n) )

    @constraint( bospp, cte[i=1:inst.m], sum(x[j] for j in inst.Cover[i]) <= 1)

    println("Solving...")

    infos = vSolve( bospp, method=:bc, verbose=false )
    println(infos)

    # ---- Querying the results
    println("Querying...")
    Y_N = getY_N( bospp )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( bospp )
    println("length X_E = ", length(X_E))

    # ---- Writing the results
    writeResults(inst.n, inst.m, fname, outputName, :bc, Y_N, X_E; infos)
end

solveDichotomy(ARGS[1])
# solve_epsilon(ARGS[1])
# solveBOBB(ARGS[1])
# solveBOBC(ARGS[1])