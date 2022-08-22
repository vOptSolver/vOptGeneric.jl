# Bi-objective set partitionning problem (2-SPA)


# ---- Packages to use
using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("../../../src/BO01BB/displayGraphic.jl")
using .vOptGeneric


function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N, X_E; total_time=nothing, infos=nothing)

    fout = open(outputName, "w")
    println(fout, "vars = $vars ; constr = $constr ")
  
    if method == :bb || method == :bc
        println(fout, infos)
    else
      println(fout, "total_times_used = $total_time")
    end
    println(fout, "size_Y_N = ", length(Y_N))
    println(fout, "Y_N = ", Y_N)
    println(fout)
    println(fout, "size_X_E = ", length(X_E))
  
    close(fout)
  
    displayGraphics(fname,Y_N, outputName)
end


# ---- Parser reading an instance of 2-SPA (format of instances compliant with vOptLib)
function loadInstance2SPA(fname::String)
    println("reading ", fname)
    f = open(fname)
    nbctr, nbvar = parse.(Int, split(readline(f))) # number of constraints , number of variables
    A = zeros(Int, nbctr, nbvar)                   # matrices of constraints
    c1 = zeros(Int, nbvar)                         # 1st vector of costs
    c2 = zeros(Int, nbvar)                         # 2nd vector of costs
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c1[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                c2[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                A[j,i] = 1
            end
        end
    end
    close(f)
    return c1, c2, A
end


# ---- Compute the set of non-dominated points Y_N of a 2SPA with vOptGeneric
function computeYNfor2SPA(  nbvar::Int,
                            nbctr::Int,
                            A::Array{Int,2},
                            c1::Array{Int,1},
                            c2::Array{Int,1},
                            method, fname; step=0.5
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


    (method == :bb || method == :bc) ? writeResults(nbvar, nbctr, "setPartitionning", fname, method, Y_N, X_E; infos) : 
                    writeResults(nbvar, nbctr, "setPartitionning", fname, method, Y_N, X_E; total_time)

end


# ---- Entry point
function main(fname::String)

    print("Compute Y_N with vOPtGeneric for a set partitionning problem with 2 objectives (2-SPA) \n\n")

    # load a numerical instance of 2SPA ----------------------------------------
    c1, c2, A = loadInstance2SPA(fname)
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    folder = "../../results/smallExamples"
    for method in [:bb, :bc] # :dicho, 
        result_dir = method≠:bb ? folder * "/" * string(method) : folder * "/" * string(method) * "/default"
            if !isdir(result_dir)
                    mkdir(result_dir)
            end
            fname = result_dir * "/" * "setPartitionning"

            computeYNfor2SPA(nbvar, nbctr, A, c1, c2, method, fname)
    end  
    
    # for step in ["0.1", "0.5", "1"]
    #     run_epsilon_ctr(step)
    # end
end


function run_epsilon_ctr(epsilon::String)
    fname = "biodidactic5.txt"

    print("Compute Y_N with vOPtGeneric for a set partitionning problem with 2 objectives (2-SPA) \n\n")

    # load a numerical instance of 2SPA ----------------------------------------
    c1, c2, A = loadInstance2SPA(fname)
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    step = parse(Float64, epsilon)
    folder = "../../results/smallExamples"
    method = :epsilon
    result_dir = folder * "/" * string(method) * "/" * string(method) * "_" * string(step) * "/"
    if !isdir(result_dir)
            mkdir(result_dir)
    end
    fname = result_dir * "setPartitionning"
    computeYNfor2SPA(nbvar, nbctr, A, c1, c2, method, fname; step=step)
end

# run_epsilon_ctr(ARGS[1])


# ---- Example
main("biodidactic5.txt")