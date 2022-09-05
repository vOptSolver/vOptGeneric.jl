

include("../../../src/BO01BB/displayGraphic.jl")


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