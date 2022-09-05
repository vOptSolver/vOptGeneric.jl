
include("../../../src/BO01BB/displayGraphic.jl")

mutable struct BOSCP
    n::Int64
    m::Int64
    C1::Vector{Int64}
    C2::Vector{Int64}
    Cover::Vector{Vector{Int64}}
end

function BOSCP()
    return BOSCP(0, 0, Vector{Int64}(), Vector{Int64}(), Vector{Vector{Int64}}())
end

function readingBOSCP(fname::String)
    f=open(fname)
    lines = readlines(f)
    close(f)

    if length(lines) == 1
        split_line = split(lines[1], "\r")
    else
        split_line = lines
    end
    # println(split_line)

    # m, n 
    inst = BOSCP()
    line = split(popfirst!(split_line), " ")
    for num in line 
        if String(num) != ""
            if inst.m < 1 
                inst.m = parse(Int64, String(num)) ; continue
            end
            if inst.n < 1 
                inst.n = parse(Int64, String(num)) ; continue
            end
        end
    end 

    # C1, C2
    for _ in 1:(inst.n/10)
        line = split(popfirst!(split_line), " ") 
        for ss in line
            if ss != ""
                push!(inst.C1, parse(Int64, String(ss)) )
            end
        end
    end

    for _ in 1:(inst.n/10)
        line = split(popfirst!(split_line), " ") 
        for ss in line
            if ss != ""
                push!(inst.C2, parse(Int64, String(ss)) )
            end
        end
    end
    @assert length(inst.C1) == inst.n ; @assert length(inst.C2) == inst.n

    # matrix 
    for i in 1:inst.m
        num = parse(Int64, split(popfirst!(split_line), " ")[end] ) ; line = split(popfirst!(split_line), " ")
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
  
    close(fout)
  
    displayGraphics(fname,Y_N, outputName)
end
