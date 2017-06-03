# -------------------------------------------------------------
# Step 1:  select JuMOO and the IP solver to use

using JuMOO
using GLPK, GLPKMathProgInterface

# -------------------------------------------------------------
# Step 2:  define the 2-objectives IP to solve

# --- Indexes, data, variables ---

j  = 1,...,5
p1 = [ 6, 4, 4, 4, 3]
p2 = [ 12, 10, 5, 3, 1]
w  = [ 8, 6, 4, 3, 2 ]

# -- Set the model ---

mKnapsack = MultiModel(solver = GLPKSolverMIP())
@variable(mKnapsack, x[1:5], Bin)
@addobjective(mKnapsack, Max , sum( p1[j] * x[j] for j=1:5 ) )
@addobjective(mKnapsack, Max , dot( p2 , x ) )
@constraint(mKnapsack, sum( w[j] * x[j] for j = 1:5 ) <= C )

# -------------------------------------------------------------
# Step 3:  call the solver (e-constraint) with the parameters

status = solve(mKnapsack, method=:eps, step=0.1) 