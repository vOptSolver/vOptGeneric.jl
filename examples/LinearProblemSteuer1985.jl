# Bi-objective linear problem (bilp)
#
# Example 6.3 (from Steuer, 1985), page 154 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using vOptGeneric, JuMP, GLPK


# ---- setting the model + values
bilp = vModel( GLPK.Optimizer )

@variable( bilp, x1 >= 0 )
@variable( bilp, x2 >= 0 )

@addobjective( bilp, Min, 3*x1 + x2 )
@addobjective( bilp, Min, -x1 - 2*x2 )

@constraint( bilp, cst1, x2 <= 3 )
@constraint( bilp, cst2, 3*x1 - x2 <= 6 )


# ---- Invoking the solver (lexicographic method)
vSolve( bilp, method=:lex )


# ---- Querying the results
Y_N = getY_N( bilp )


# ---- Displaying the results for lex(1,2) and lex(2,1)
for s = 1:2
    print("X = [")
    print(round(getvOptData(bilp).X_E[s][1],digits=5))
    print(" ")
    print(round(getvOptData(bilp).X_E[s][2],digits=5))
    print("] | Z = [")
    print(round(Y_N[s][1], digits=5))
    print(" ")
    print(round(Y_N[s][2], digits=5))
    println("]")
end
