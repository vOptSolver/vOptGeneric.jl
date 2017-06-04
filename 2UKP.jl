# =============================================================
# vOptSolver - June 2017
#
# Example where:
# - the problem is a bi-objective 0/1 unidimensional knapsack problem
# - the data are provided explicitely using the Julia syntax
# - the problem is defined with the algebraic language (=> it is viewed as a non-structured problem)
# - the algorithm is the epsilon-contraint method

# -------------------------------------------------------------
# Step 1:  announce the use of (1) JuMOO and (2) the IP solver to call

using JuMOO
using GLPK, GLPKMathProgInterface

# -------------------------------------------------------------
# Step 2:  define the 2-objectives IP (here: 2UKP) to solve

# --- Index, data ---

j  = 1,...,5            # index on the 5 items to consider
p1 = [ 6, 4, 4, 4, 3]   # profits of items for the objective 1
p2 = [ 12, 10, 5, 3, 1] # profits of items for the objective 2
w  = [ 8, 6, 4, 3, 2 ]  # weight of items for the constraint
C  = 12                 # capacity of the knapsack

# -- Set the model (with an explicit and an implicit formulation of the objectives) ---

mKnapsack = MultiModel(solver = GLPKSolverMIP())
@variable(mKnapsack, x[1:5], Bin)
@addobjective(mKnapsack, Max , sum( p1[j] * x[j] for j=1:5 ) )
@addobjective(mKnapsack, Max , dot( p2 , x ) )
@constraint(mKnapsack, sum( w[j] * x[j] for j = 1:5 ) <= C )

# -------------------------------------------------------------
# Step 3:  call the solver (here: e-constraint) with the parameters

status = solve(mKnapsack, method=:eps, step=0.1) 

# -------------------------------------------------------------
# Step 4:  Get the results
