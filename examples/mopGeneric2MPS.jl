# Example of vOptGeneric, for a 2IP (the 2ap03.mop is a 2LAP)
# the numerical instance is read on a file, according to the MOP format
# July 2017-updated July 2021

# ---- Packages to use
using vOptGeneric, JuMP, GLPK

m = parseMOP( "2ap03.mop" , GLPK.Optimizer ) # comment: set up the path to the file if required

# ---- Invoking the solver (Chalmet method)
@time vSolve( m, method=:chalmet, step=0.5 )

printX_E(m)
getY_N(m)
