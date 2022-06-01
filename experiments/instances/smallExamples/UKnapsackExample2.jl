# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Example 2 with 20 variables


# ---- Packages to use
using JuMP, GLPK, LinearAlgebra

include("../../../src/vOptGeneric.jl")
include("../../../src/BO01BB/displayGraphic.jl")
using .vOptGeneric


function writeResults(fname::String, outputName::String, method, Y_N; total_time=nothing, infos=nothing)

        fout = open(outputName, "w")
      
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


function BOUKP(method, fname)

    # ---- Values of the instance to solve
    p1 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62]
    p2 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
    w  = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99]
    c  = 900
    size = length(p1)


    # ---- setting the model
    m = vModel( GLPK.Optimizer ) ; JuMP.set_silent(m)
    @variable( m, x[1:size], Bin )
    @addobjective( m, Max, dot(x, p1) )
    @addobjective( m, Max, dot(x, p2) )
    @constraint( m, dot(x, w) <= c )


    if method == :bb
        infos = vSolve( m, method=:bb, verbose=false )
    elseif method == :dicho 
        start = time()
        vSolve( m, method=:dicho, verbose=false )
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( m, method=:epsilon, step=0.5, verbose=false )
        total_time = round(time() - start, digits = 2)
    else
        @error "Unknown method parameter $(method) !"
    end

    # ---- Querying the results
    Y_N = getY_N( m )

    (method == :bb) ? writeResults("UKnapsackExample2", fname, method, Y_N; infos) : 
                    writeResults("UKnapsackExample2", fname, method, Y_N; total_time)

end


function main()
    folder = "../../results/smallExamples/"
    for method in [:bb, :dicho, :epsilon]
            result_dir = folder * "/" * string(method)
            if !isdir(result_dir)
                    mkdir(result_dir)
            end
            fname = result_dir * "/" * "UKnapsackExample2"

            # # if method == :bb 
                BOUKP(method, fname) 
            # end
    end
end


main()
