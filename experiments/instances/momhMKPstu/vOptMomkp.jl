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
using JuMP, GLPK
include("../../../src/vOptGeneric.jl")
using .vOptGeneric

# ==============================================================================
# Setting a Bi01IP instance and solving with vOptGeneric
function vSolveBi01IP(solverSelected, C, A, B)

  m, n = size(A)

  # ---- setting the model
  println("Building...")
  Bi01IP = vModel( solverSelected )
  @variable( Bi01IP, x[1:n], Bin )
  @addobjective( Bi01IP, Max, sum(C[1,j] * x[j] for j=1:n) )
  @addobjective( Bi01IP, Max, sum(C[2,j] * x[j] for j=1:n) )
  @constraint( Bi01IP, cte[i=1:m], sum(A[i,j] * x[j] for j=1:n) <= B[i])

  # ---- Invoking the solver (epsilon constraint method)
  println("Solving...")
  # vSolve( Bi01IP, method=:dicho, verbose=false )
  # vSolve( Bi01IP, method=:epsilon, step=0.5, verbose=false )
  infos = vSolve( Bi01IP, method=:bb, verbose=true )
  println(infos)



  # ---- Querying the results
  println("Querying...")
  Y_N = getY_N( Bi01IP )

  return Y_N
end

# ==============================================================================

# # Example on how to use it :

# solverSelected = GLPK.Optimizer
# Y_N, X_E = setBi01IP(solverSelected, C, A, B)

# # ---- Displaying the results (X_E and Y_N)
# for n = 1:length(Y_N)
#    X = value.(X_E, n)
#    print(findall(elt -> elt ≈ 1, X))
#    println("| z = ",Y_N[n])
# end

# for n = 1:length(Y_N)
#    println(Y_N[n])
# end
