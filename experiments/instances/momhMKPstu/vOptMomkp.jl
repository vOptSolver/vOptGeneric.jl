# ==============================================================================
# Xavier Gandibleux - November 2021
#   Implemented in Julia 1.6

# ==============================================================================
# Using ϵ-constraint / dichotomy / branch-and-bound method with vOptGeneric, compute the set of non-dominated
# points of the following problem:
#
#   Max  sum{j=1,...,n} C(k,j) x(j)                k=1,...,2
#   st   sum{j=1,...,n} A(i,j) x(j) <= B(i)        i=1,...,m
#                       x(j) = 0 or 1

# ==============================================================================

println("""\nProject MOMH 2021" --------------------------------\n""")

const verbose = false
const graphic = true

include("parserMomkpZL.jl")
include("parserMomkpPG.jl")

using JuMP, CPLEX
include("../../../src/vOptGeneric.jl")
include("../../../src/BO01BB/displayGraphic.jl")
using .vOptGeneric


# ==============================================================================
# Datastructure of a generic bi-objective 0/1 IP where all coefficients are integer
#
#   Max  sum{j=1,...,n} C(k,j) x(j)                k=1,...,p
#   st   sum{j=1,...,n} A(i,j) x(j) <= b_{i}       i=1,...,m
#                       x(j) = 0 or 1

struct _bi01IP
  C  :: Matrix{Int} # objective functions, k=1..2, j=1..n
  A  :: Matrix{Int} # matrix of constraints, i=1..m, j=1..n
  b  :: Vector{Int} # right-hand side, i=1..m
end


function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N, X_E; total_time=nothing, infos=nothing)

  fout = open(outputName, "w")
  println(fout, "vars = $vars ; constr = $constr ")

  if method == :bb
    println(fout, infos)
  else
    println(fout, "total_times_used = $total_time")
  end
  println(fout, "size_Y_N = ", length(Y_N))
  println(fout, "Y_N = ", Y_N)
  println(fout)
  println(fout, "size_X_E = ", length(X_E))
  # println(fout, "X_E = ", X_E)

  close(fout)

  displayGraphics(fname,Y_N, outputName)
end


function vSolveBi01IP(solverSelected, C, A, B, fname, method)

  println("method : ", string(method))

  result_dir = "../../results/momhMKPstu/" * split(fname, "/")[2]
  if !isdir(result_dir)
    mkdir(result_dir)
  end

  folder = result_dir * "/" * string(method)
  if !isdir(folder)
    mkdir(folder)
  end

  m, n_before = size(A)
  # scale test
  for n = 10:10:50

    ratio = n/n_before

    subfolder = folder * "/" * string(n)
    if !isdir(subfolder)
      mkdir(subfolder)
    end

    outputName = subfolder * "/" * split(fname, "/")[end]
    # # TODO : if a file already exists
    # if isfile(outputName) # && method != :bb 
    #   return
    # end

    # ---- setting the model
    println("Building...")
    Bi01IP = vModel( solverSelected ) ; JuMP.set_silent(Bi01IP)
    @variable( Bi01IP, x[1:n], Bin )
    @addobjective( Bi01IP, Max, sum(C[1,j] * x[j] for j=1:n) )
    @addobjective( Bi01IP, Max, sum(C[2,j] * x[j] for j=1:n) )
    @constraint( Bi01IP, cte[i=1:m], sum(A[i,j] * x[j] for j=1:n) <= round(Int, B[i]*ratio))

    # ---- Invoking the solver (epsilon constraint method)
    println("Solving...")
    if method == :dicho
      start = time()
      vSolve( Bi01IP, method=:dicho, verbose=false)
      total_time = round(time() - start, digits = 2)
    elseif method == :epsilon
      start = time()
      vSolve( Bi01IP, method=:epsilon, step=0.5, verbose=false )
      total_time = round(time() - start, digits = 2)
    elseif method == :bb
      infos = vSolve( Bi01IP, method=:bb, verbose=false )
      println(infos)
    end

    # ---- Querying the results
    println("Querying...")
    Y_N = getY_N( Bi01IP )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( Bi01IP )
    println("length X_E = ", length(X_E))

    # ---- Writing the results
    (method == :bb) ?
      writeResults(n, m, fname, outputName, method, Y_N, X_E; infos) : 
      writeResults(n, m, fname, outputName, method, Y_N, X_E; total_time)

    # return Y_N
  end

end


function main(fname::String)
  if !isfile(fname)
    @error "This file doesn't exist ! $fname"
  end

  if split(fname, "/")[end] == "readme" return end

  inst = split(fname, "/")[2]
  if inst == "MOMKP"
    # Read and load an instance of MO-01MKP from the collection of E. Zitzler / M. Laumanns
    momkpZL = readInstanceMOMKPformatZL(verbose,fname)

    # Reduce the MO-MKP instance to the two first objectives and store it as a generic Bi-01IP
    dat = _bi01IP(momkpZL.P[1:2,:], momkpZL.W, momkpZL.ω)
  elseif inst == "MOBKP"
    # Read and load an instance of Bi-01BKP from the collection of O. Perederieieva / X. Gandibleux
    momkpPG = readInstanceMOMKPformatPG(verbose,fname)

    # Store it as a generic bi-01IP
    dat = _bi01IP(momkpPG.P, momkpPG.W, momkpPG.ω)
  else
    @error "Unknown input file $fname"
  end

  solverSelected = CPLEX.Optimizer
  for method in [:dicho, :epsilon, :bb] # 
    vSolveBi01IP(solverSelected, dat.C, dat.A, dat.b, fname, method) 
  end

end


main(ARGS[1])
