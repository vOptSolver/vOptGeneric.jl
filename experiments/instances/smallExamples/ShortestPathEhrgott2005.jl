# Bi-objective shortest path problem (bisp)
#
# Exercise 9.5 page 269 of
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

function vSolveBOSP(method, fname; step=0.5)
        # ---- Values of the instance to solve
        M  = 50
        C1 =  [ M 4 5 M M M ;  # coefficients's vector for arc (i,j) of the objective 1
                M M 2 1 2 7 ;
                M M M 5 2 M ;
                M M 5 M M 3 ;
                M M M M M 4 ;
                M M M M M M  ]

        C2 =  [ M 3 1 M M M ; # coefficients's vector for arc (i,j) of the objective 2
                M M 1 4 2 2 ;
                M M M 1 7 M ;
                M M 1 M M 2 ;
                M M M M M 2 ;
                M M M M M M  ]

        n  = size(C2,1)       # number of nodes

        # ---- setting the model
        bisp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent( bisp )

        @variable( bisp, x[1:n,1:n], Bin )

        @addobjective( bisp, Min , sum(x[i,j]*C1[i,j] for i=1:n, j=1:n) )
        @addobjective( bisp, Min , sum(x[i,j]*C2[i,j] for i=1:n, j=1:n) )

        @constraint( bisp, node_s, sum(x[1,j] for j = 1:n) == 1 )
        @constraint( bisp, node_t, sum(x[i,n] for i = 1:n) == 1 )
        @constraint( bisp, cstr[i=2:n-1], sum(x[i,j] for j = 1:n) - sum(x[j,i] for j =1:n) == 0 )


        if method == :bb
                infos = vSolve( bisp, method=:bb, verbose=false )
        elseif method == :dicho 
                start = time()
                vSolve( bisp, method=:dicho, verbose=false )
                total_time = round(time() - start, digits = 2)
        elseif method==:epsilon 
                start = time()
                vSolve( bisp, method=:epsilon, step=step, verbose=false )
                total_time = round(time() - start, digits = 2)
        else
                @error "Unknown method parameter $(method) !"
        end

        # ---- Querying the results
        Y_N = getY_N( bisp )

        (method == :bb) ? 
                writeResults(n*n, n, "ShortestPathEhrgott2005", fname, method, Y_N; infos) : 
                writeResults(n*n, n, "ShortestPathEhrgott2005", fname, method, Y_N; total_time)

end


function main()
        folder = "../../results/smallExamples/"
        for method in [:dicho, :bb] # 
                result_dir = method≠:bb ? folder * "/" * string(method) : folder * "/" * string(method) * "/default"
                if !isdir(result_dir)
                        mkdir(result_dir)
                end
                fname = result_dir * "/" * "ShortestPathEhrgott2005"

                vSolveBOSP(method, fname)
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
        fname = result_dir * "/" * "ShortestPathEhrgott2005"

        vSolveBOSP(method, fname; step=step)
end

# run_epsilon_ctr(ARGS[1])

main()
