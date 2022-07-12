# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Exercise 10.2 page 290 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("../../../src/BO01BB/displayGraphic.jl")
using .vOptGeneric



function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N; total_time=nothing, infos=nothing)

    fout = open(outputName, "w")
    println(fout, "vars = $vars ; constr = $constr ")
  
    if method == :bb
      println(fout, infos)
    else
      println(fout, "total_times_used = $total_time")
    end
    println(fout, "size_Y_N = ", length(Y_N))
    println(fout, "Y_N = ", Y_N)
  
    close(fout)
  
    displayGraphics(fname,Y_N, outputName)
end

function vSolveBOUKP(method, fname; step=0.5)
    
    # ---- Values of the instance to solve
    p1 = [10, 3,  6, 8, 2]  # coefficients's vector of the objective 1
    p2 = [12, 9, 11, 5, 6]  # coefficients's vector of the objective 2
    w  = [ 4, 5,  2, 5, 6]  # coefficients's vector of weights
    c  = 17                 # nominal capacity
    n  = length( p1 )       # number of items


    # ---- setting the model
    biukp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent( biukp )
    @variable( biukp, x[1:n], Bin )
    @addobjective( biukp, Max, sum( p1[j]*x[j] for j=1:n ) )
    @addobjective( biukp, Max, sum( p2[j]*x[j] for j=1:n ) )
    @constraint( biukp, sum( w[j]*x[j] for j=1:n ) <= c )

    if method == :bb
        infos = vSolve( biukp, method=:bb, verbose=false )
    elseif method == :dicho 
        start = time()
        vSolve( biukp, method=:dicho, verbose=false )
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( biukp, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)
    else
        @error "Unknown method parameter $(method) !"
    end

    # ---- Querying the results
    Y_N = getY_N( biukp )

    (method == :bb) ? 
        writeResults(n, 1, "UKnapsackEhrgott2005", fname, method, Y_N; infos) : 
        writeResults(n, 1, "UKnapsackEhrgott2005", fname, method, Y_N; total_time)
end


function main()
    folder = "../../results/smallExamples/"
    for method in [:dicho, :bb] # 
        result_dir = method≠:bb ? folder * "/" * string(method) : folder * "/" * string(method) * "/default"
            if !isdir(result_dir)
                    mkdir(result_dir)
            end
            fname = result_dir * "/" * "UKnapsackEhrgott2005"

            vSolveBOUKP(method, fname)
    end

    for step in ["0.1", "1", "5"]
        run_epsilon_ctr(step)
    end
end


function run_epsilon_ctr(epsilon::String)
    step = parse(Float64, epsilon)
    folder = "../../results/smallExamples/"
    method = :epsilon
    result_dir = folder * "/" * string(method) * "/" * string(method) * "_" * string(step)
    if !isdir(result_dir)
            mkdir(result_dir)
    end
    fname = result_dir * "/" * "UKnapsackEhrgott2005"

    vSolveBOUKP(method, fname; step=step)  
end

# run_epsilon_ctr(ARGS[1])

main()