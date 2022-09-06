# ---- Packages to use
using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("parserBOSPA.jl")
using .vOptGeneric


# ---- Compute the set of non-dominated points Y_N of a 2SPA with vOptGeneric
function computeYNfor2SPA(  nbvar::Int,
    nbctr::Int,
    A::Array{Int,2},
    c1::Array{Int,1},
    c2::Array{Int,1},
    method, fname, outputName; step=0.5
 )

    # ---- setting the model
    model = vModel( CPLEX.Optimizer )
    JuMP.set_silent(model)

    @variable(model, x[1:nbvar], Bin)
    @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    @addobjective(model, Min, sum(c1[i]*x[i] for i in 1:nbvar))
    @addobjective(model, Min, sum(c2[i]*x[i] for i in 1:nbvar))

    if method == :bb
        infos = vSolve( model, method=:bb, verbose=false )
    elseif method == :bc 
        infos = vSolve( model, method=:bc, verbose=false )
    elseif method == :dicho 
        start = time()
        vSolve( model, method=:dicho, verbose=false )
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( model, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)
    else
        @error "Unknown method parameter $(method) !"
    end

    # ---- Querying the results
    Y_N = getY_N( model )

    X_E = getX_E( model )


    (method == :bb || method == :bc) ? writeResults(nbvar, nbctr, fname, outputName, method, Y_N, X_E; infos) :
        writeResults(nbvar, nbctr, fname, outputName, method, Y_N, X_E; total_time)

end


function solve(fname::String, method)

    # load a numerical instance of 2SPA ----------------------------------------
    c1, c2, A = loadInstance2SPA(fname)
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2
    if nbvar >= 1000 return end

    folder = "../../results/SPA/BOSPA"
    if !isdir(folder)
        mkdir(folder)
    end
    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end
    inst_name = split(fname, "/")[end]

    println("n=$nbvar m=$nbctr ") ; outputName = result_folder * "/" * split(inst_name, ".")[1] * ".dat"
    if isfile(outputName) return end #TODO : ignore existed file  
    computeYNfor2SPA(nbvar, nbctr, A, c1, c2, method, string(split(inst_name, ".")[1]), outputName)
end

solve(ARGS[1], :dicho)
solve(ARGS[1], :epsilon)
solve(ARGS[1], :bb)
solve(ARGS[1], :bc)
