# Bi-objective linear assignment problem (bilap)
#
# Example 9.38 (from Ulungu and Teghem, 1994), page 255 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using JuMP, GLPK

include("../../src/vOptGeneric.jl")
using .vOptGeneric


# ---- Values of the instance to solve
C1 = [  5  1  4  7 ;  # coefficients's vector of the objective 1
        6  2  2  6 ;
        2  8  4  4 ;
        3  5  7  1   ]

C2 = [  3  6  4  2 ;  # coefficients's vector of the objective 2
        1  3  8  3 ;
        5  2  2  3 ;
        4  2  3  5   ]

n  = size(C2,1)       # number of lines/columns


# ---- setting the model
bilap = vModel( GLPK.Optimizer )

@variable( bilap, x[1:n,1:n], Bin )

@addobjective( bilap , Min, sum( C1[i,j]*x[i,j] for i=1:n,j=1:n ) )
@addobjective( bilap , Min, sum( C2[i,j]*x[i,j] for i=1:n,j=1:n ) )

@constraint( bilap , cols[i=1:n], sum(x[i,j] for j=1:n) == 1 )
@constraint( bilap , rows[j=1:n], sum(x[i,j] for i=1:n) == 1 )


# ---- Invoking the solver (branch and bound method)
vSolve( bilap, method=:bb, verbose=true )


# # ---- Invoking the solver (Chalmet method)
# vSolve( bilap, method=:chalmet, step=0.5, verbose=true )

# # ---- Querying the results
# Y_N = getY_N( bilap )

# # ---- Displaying the results
# printX_E( bilap )



# method=:chalmet, Output : 

# [6.0, 24.0] : x[3,1]=1 x[1,2]=1 x[2,3]=1 x[4,4]=1 
# [9.0, 17.0] : x[3,1]=1 x[2,2]=1 x[1,3]=1 x[4,4]=1 
# [12.0, 13.0] : x[1,1]=1 x[2,2]=1 x[3,3]=1 x[4,4]=1 
# [16.0, 11.0] : x[4,1]=1 x[2,2]=1 x[3,3]=1 x[1,4]=1 
# [19.0, 10.0] : x[2,1]=1 x[4,2]=1 x[1,3]=1 x[3,4]=1 
# [22.0, 7.0] : x[2,1]=1 x[4,2]=1 x[3,3]=1 x[1,4]=1 

# method=:bb, Output :

# incumbent : NaturalOrderVector[
# Solution( 
#  xEquiv = [[0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]]
#  y = [22.0, 7.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
#  y = [12.0, 13.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
#  y = [9.0, 17.0]
#  is_binary ? true )
# Solution( 
#  xEquiv = [[0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
#  y = [6.0, 24.0]
#  is_binary ? true )
# ]
