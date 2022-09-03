# =============================================================================
"""
    0/1 uncapacitated facility location problem with two objectives (biUFLP)
    uncapacitatedFacilityLocation.jl
    September 2022
"""


# =============================================================================

print("Loading and compiling vOptGeneric, JuMP, GLPK...")
using vOptGeneric, JuMP, GLPK, Printf
println(" done!")


# =============================================================================
# structure of a 2UFLP instance

mutable struct Instance
    fname :: String        # name of the file
    nI    :: Int64         # number of users
    nJ    :: Int64         # number of services
    c1    :: Array{Int,2}  # assignment costs users-services for objective 1
    c2    :: Array{Int,2}  # assignment costs users-services for objective 2
    r1    :: Array{Int,1}  # running costs for objective 1
    r2    :: Array{Int,1}  # running costs for objective 2
end


# =============================================================================
"""
    load2UFLP(fdirectory::String, fname::String)

Load an instance of bi-objective UFLP
-  ↓  fdirectory : path to the data file
-  ↓  fname      : file name of the data
-  ↑             : data of an instance
"""

function load2UFLP(fdirectory::String, fname::String)

    f=open(fdirectory*"/"*fname)
 
    # read the number of users (nI)
    nI = parse(Int, readline(f) )
    # read the number of services (nJ) 
    nJ = parse(Int, readline(f) )
    # read the following line (separator -> line without information)
    useless = readline(f)

    # assignment costs users-services (matrix nI x nJ) ------------------------
    c1 = Array{Int64,2}(undef,nI,nJ)
    c2 = Array{Int64,2}(undef,nI,nJ)

    # objective 1
    for i=1:nI
        c1[i,:] = parse.(Int64, split(readline(f)) )
    end
    # read the following line (separator -> line without information)
    useless = readline(f)

    # objective 2
    for i=1:nI
        c2[i,:] = parse.(Int64, split(readline(f)) )
    end
    # read the following line (separator -> line without information)
    useless = readline(f)

    # running costs of services (vector nJ) -----------------------------------
    r1 = Array{Int64,1}(undef,nJ)
    r2 = Array{Int64,1}(undef,nJ)

    # objective 1
    r1[:] = parse.(Int64, split(readline(f)) )
    # read the following line (separator -> line without information)
    useless = readline(f)

    # objective 2
    r2[:] = parse.(Int64, split(readline(f)) )
   
    close(f)

    return Instance(fname,nI,nJ,c1,c2,r1,r2)

end


# ==============================================================================
"""
    createModel2UFLP(solver::DataType, data::Instance)

    Create the vOptGeneric model of 2UFLP 
"""
function createModel2UFLP(solver::DataType, data::Instance)

    model = vModel( solver )
    
    @variable(model, x[1:data.nI,1:data.nJ], Bin)
    @variable(model, s[1:data.nJ], Bin)

    @addobjective( model, Min, sum(data.c1[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r1[j]*s[j] for j in 1:data.nJ) )
    @addobjective( model, Min, sum(data.c2[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r2[j]*s[j] for j in 1:data.nJ) )

    @constraint( model, [i=1:data.nI], sum(x[i,j] for j in 1:data.nJ) == 1 )
    @constraint( model, [i=1:data.nI,j=1:data.nJ], x[i,j] <= s[j] )    

    return model
end


# ==============================================================================
function main()

    # -------------------------------------------------------------------------
    # Load an instance (files are available in vOptLib)

    dname = "."               # path to the folder containing the files
    fname = "didactic1.txt"   # filename of the instance to solve
    data  = load2UFLP(dname,fname)

    # -------------------------------------------------------------------------
    # Display information about the instance to solve

    println("\nfilename      : $(data.fname)") 
    println("nI (users)    : $(data.nI)") 
    println("nJ (services) : $(data.nJ)\n") 

    # -------------------------------------------------------------------------
    # compute YN with vOptGeneric using ϵ-constraint method and GLPK

    println("Running the ϵ-constraint method with GLPK... \n") 
    solver = GLPK.Optimizer    
    mod2UFLP = createModel2UFLP(solver, data)
    set_silent(mod2UFLP)

    getTime = time()
    vSolve(mod2UFLP, method=:epsilon, step = 1.0)
    timevOPt = round(time()- getTime, digits=4)


    # -------------------------------------------------------------------------
    # Display the results 

    println("\nDisplaying the results... \n") 
    YN = getY_N( mod2UFLP )
    println("fname: $(fname[1:end-4])     tOpt: $(timevOPt) sec     #YN: $(length(YN)) points\n")
    printX_E( mod2UFLP )
    println("\n...done!")

    return nothing
end

# ==============================================================================
main()
