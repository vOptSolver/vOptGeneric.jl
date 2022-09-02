
# reading instances of Bi-objective set packing problems 

include("../../../src/BO01BB/displayGraphic.jl")

mutable struct BOSPP
    n::Int64
    m::Int64
    C1::Vector{Int64}
    C2::Vector{Int64}
    Cover::Vector{Vector{Int64}}
end

function BOSPP()
    return BOSPP(0, 0, Vector{Int64}(), Vector{Int64}(), Vector{Vector{Int64}}())
end

function readingBOSPP(fname::String)
    f=open(fname)
    lines = readlines(f)
    close(f)

    # m, n 
    inst = BOSPP()
    line = split(popfirst!(lines), " ")
    inst.m = parse(Int64, String(line[1])) ; inst.n = parse(Int64, String(line[2]))

    # C1, C2
    line = split(popfirst!(lines), " ") ; inst.C1 = [parse(Int64, String(line[i])) for i in 1:length(line) if line[i] != ""]
    line = split(popfirst!(lines), " ") ; inst.C2 = [parse(Int64, String(line[i])) for i in 1:length(line) if line[i] != ""]
    @assert length(inst.C1) == inst.n ; @assert length(inst.C2) == inst.n

    # matrix 
    for i in 1:inst.m
        num = parse(Int64, split(popfirst!(lines), " ")[end] ) ; line = split(popfirst!(lines), " ")
        push!(inst.Cover, [parse(Int64, String(line[i])) for i=1:length(line) if line[i] != ""])
    end
    return inst 
end

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
    # println(fout, "X_E = ", X_E)
  
    close(fout)
  
    displayGraphics(fname,Y_N, outputName)
end
