# Bi-objective shortest path problem (bisp)
#
# Exercise 9.5 page 269 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using vOptGeneric, JuMP, GLPK


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
bisp = vModel( GLPK.Optimizer ) ; JuMP.set_silent( bisp )

@variable( bisp, x[1:n,1:n], Bin )

@addobjective( bisp, Min , sum(x[i,j]*C1[i,j] for i=1:n, j=1:n) )
@addobjective( bisp, Min , sum(x[i,j]*C2[i,j] for i=1:n, j=1:n) )

@constraint( bisp, node_s, sum(x[1,j] for j = 1:n) == 1 )
@constraint( bisp, node_t, sum(x[i,n] for i = 1:n) == 1 )
@constraint( bisp, cstr[i=2:n-1], sum(x[i,j] for j = 1:n) - sum(x[j,i] for j =1:n) == 0 )


# ---- Invoking the solver (epsilon constraint method)
vSolve( bisp, method=:epsilon, step=0.5, verbose=true )


# ---- Querying the results
Y_N = getY_N( bisp )


# ---- Displaying the results
for s = 1:length(Y_N)
    X = value.(x, s)
    print("Path: ")
    for ind in findall(val -> val ≈ 1, X)
        i,j = ind.I
        print(" $i->$j ")
    end
    println("| Z = ",Y_N[s])
end
