# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
using Printf
varname_generic(m::Model, col::Integer) = "VAR$(col)"

function varname_given(m::Model, col::Integer)
    # TODO: deal with non-ascii characters?
    name = getname(m, col)
    for (pat, sub) in [("[", "_"), ("]", ""), (",", "_")]
        name = replace(name, pat => sub)
    end
    name
end


function writeMOP(m, fname::AbstractString, genericnames=false)
    varname = genericnames ? varname_generic : varname_given

    f = open(fname, "w")

    write(f,"NAME   vOptModel\n")

    numRows = length(m.linconstr)


    md = getvOptData(m)
    objlincoef = [prepAffObjective(m, f) for f in md.objs]

    write(f,"OBJSENSE\n")
    if md.objSenses[1] == :Min
        write(f," MIN\n")
    else
        write(f," MAX\n")
    end

    for i = 2:length(md.objs)
        if md.objSenses[i] != md.objSenses[1]
            println("Warning, MOP does not support different optimisation senses. Flipping objective coefficients.")
            objlincoef[i] = -objlincoef[i]
        end
    end

    # Objective and constraint names
    write(f,"ROWS\n")

    for i = 1:length(md.objs)
        @printf(f," N  OBJ%d\n", i)
    end

    hasrange = false
    for c in 1:numRows
        rowsense = JuMP.sense(m.linconstr[c])
        if rowsense == :(<=)
            senseChar = 'L'
        elseif rowsense == :(==)
            senseChar = 'E'
        elseif rowsense == :(>=)
            senseChar = 'G'
        else
            hasrange = true
            senseChar = 'E'
        end
        @printf(f," %c  CON%d\n",senseChar,c)
    end

    

    A = JuMP.prepConstrMatrix(m)
    colptr = A.colptr
    rowval = A.rowval
    nzval = A.nzval

    # Output each column
    inintegergroup = false
    write(f,"COLUMNS\n")
    for col in 1:m.numCols
        t = m.colCat[col]
        (t == :SemiCont || t == :SemiInt) && error("The MPS file writer does not currently support semicontinuous or semi-integer variables")
        if (t == :Bin || t == :Int) && !inintegergroup
            @printf(f,"    MARKER    'MARKER'                 'INTORG'\n")
            inintegergroup = true
        elseif (t == :Cont || t == :Fixed) && inintegergroup
            @printf(f,"    MARKER    'MARKER'                 'INTEND'\n")
            inintegergroup = false
        end
        for ind in colptr[col]:(colptr[col+1]-1)
            @printf(f,"    %s  CON%d  ",varname(m,col),rowval[ind])
            (Base.Grisu).print_shortest(f,nzval[ind])
            println(f)
        end
        for obj = 1:length(md.objs)
            @printf(f,"    %s  OBJ%d  ",varname(m,col), obj)
            (Base.Grisu).print_shortest(f,objlincoef[obj][col])
            println(f)
        end
    end
    if inintegergroup
        @printf(f,"    MARKER    'MARKER'                 'INTEND'\n")
    end

    # RHSs
    write(f,"RHS\n")
    for c in 1:numRows
        rowsense = JuMP.sense(m.linconstr[c])
        if rowsense != :range
            @printf(f,"    RHS    CON%d    ",c)
            (Base.Grisu).print_shortest(f,JuMP.rhs(m.linconstr[c]))
        else
            @printf(f,"    RHS    CON%d    ",c)
            (Base.Grisu).print_shortest(f,m.linconstr[c].lb)
        end
        println(f)
    end

    # RANGES
    if hasrange
        write(f,"RANGES\n")
        for c in 1:numRows
            rowsense = JuMP.sense(m.linconstr[c])
            if rowsense == :range
                @printf(f,"    rhs    CON%d    ",c)
                (Base.Grisu).print_shortest(f,m.linconstr[c].ub-m.linconstr[c].lb)
                println(f)
            end
        end
    end


    # BOUNDS
    write(f,"BOUNDS\n")
    for col in 1:m.numCols
        if m.colLower[col] == 0
            if m.colUpper[col] != Inf
                # Default lower 0, and an upper
                @printf(f,"  UP BOUND %s ", varname(m,col))
                (Base.Grisu).print_shortest(f, m.colUpper[col])
                println(f)
            else
                # Default bounds. Explicitly state for solvers like Gurobi. See issue #792
                @printf(f,"  PL BOUND %s", varname(m,col))
                println(f)
            end
        elseif m.colLower[col] == -Inf && m.colUpper[col] == +Inf
            # Free
            @printf(f, "  FR BOUND %s\n", varname(m,col))
        elseif m.colLower[col] != -Inf && m.colUpper[col] == +Inf
            # No upper, but a lower
            @printf(f, "  PL BOUND %s\n  LO BOUND %s ",varname(m,col),varname(m,col))
            (Base.Grisu).print_shortest(f,m.colLower[col])
            println(f)
        elseif m.colLower[col] == -Inf && m.colUpper[col] != +Inf
            # No lower, but a upper
            @printf(f,"  MI BOUND %s\n  UP BOUND %s ",varname(m,col),varname(m,col))
            (Base.Grisu).print_shortest(f,m.colUpper[col])
            println(f)
        else
            # Lower and upper
            @printf(f, "  LO BOUND %s ",varname(m,col))
            (Base.Grisu).print_shortest(f,m.colLower[col])
            println(f)
            @printf(f, "  UP BOUND %s ",varname(m,col))
            (Base.Grisu).print_shortest(f,m.colUpper[col])
            println(f)
        end
    end
    
    write(f,"ENDATA\n")
    close(f)
    nothing
end

nextline(f) = split(chomp(readline(f)), ' ', keepempty=false)
function parseMOP(fname::AbstractString; solver=JuMP.UnsetSolver())

    m = vModel(solver = solver)

    open(fname) do f
        ln = nextline(f)
        while ln[1] != "NAME"
            ln = nextline(f)
        end

        objSense = :Min
        ln = nextline(f)
        if ln[1] == "OBJSENSE"
            ln = nextline(f)
            if ln[1] == "MAXIMIZE" || ln[1] == "MAX"
                objSense = :Max
            end
            ln = nextline(f)
        end

        if ln[1] == "ROWS"
            ln = nextline(f)
            DicoExpr = Dict{AbstractString,Tuple{JuMP.AffExpr, Symbol}}()
            objOrder = Vector{AbstractString}()
            while ln[1] != "COLUMNS"
                push!(DicoExpr, ln[2] => ( zero(JuMP.AffExpr), Symbol(ln[1]) ) )
                ln[1] == "N" && push!(objOrder, ln[2])
                ln = nextline(f)
            end
        end

        if ln[1] == "COLUMNS"
            ln = nextline(f)
            DicoVar = Dict{AbstractString,JuMP.Variable}()
            isInt = false
            while ln[1] != "RHS"

                if ln[3] == "'INTORG'"
                    isInt = true
                elseif ln[3] == "'INTEND'"
                    isInt = false
                else
                    if !(ln[1] in keys(DicoVar))
                        
                        _cat = isInt ? (:Int) : (:Cont)
                        v = @variable(m, category=_cat, basename=ln[1], lowerbound=0.0)

                        # v = @variable(m, basename=ln[1], lowerbound=0.0)
                        # if isInt 
                        #     m.colCat[v.col] = :Int
                        # end

                        push!(DicoVar, ln[1] => v)
                    end

                    expr, expr_type = DicoExpr[ln[2]]
                    var = DicoVar[ln[1]]
                    coeff = parse(Float64, ln[3])

                    push!(expr, coeff, var)

                end
                ln = nextline(f)
            end
        end

        DicoCstr = Dict{AbstractString, Tuple{JuMP.LinearConstraint, Symbol}}()
        DicoObj  = Dict{AbstractString, JuMP.QuadExpr}()
        for (k,v) in DicoExpr
            expr, expr_type = v
            if expr_type == :N
                push!(DicoObj, k => convert(QuadExpr, expr))
            else
                push!(DicoCstr, k => (LinearConstraint(expr, -Inf, Inf), expr_type) )
            end
        end

        md = m.ext[:vOpt]
        for key in objOrder
            obj = DicoObj[key]
            push!(md.objs, obj)
            push!(md.objSenses, objSense)
        end

        if ln[1] == "RHS"
            ln = nextline(f)
            while length(ln) > 1
                
                constraint, constraint_type = DicoCstr[ln[2]]
                if constraint_type == :E
                    constraint.lb = parse(Float64, ln[3])
                    constraint.ub = parse(Float64, ln[3])
                elseif constraint_type == :G
                    constraint.lb = parse(Float64, ln[3])
                elseif constraint_type == :L
                    constraint.ub = parse(Float64, ln[3])
                end
                ln = nextline(f)
            end
        end

        for (k,v) in DicoCstr
            cstr, cstr_type = v
            JuMP.addconstraint(m, cstr)
        end


        if ln[1] == "RANGES"
            ln = nextline(f)
            while length(ln) > 1
                
                constraint, constraint_type = DicoCstr[ln[2]]
                r = parse(Float64, ln[3])
                if constraint_type == :G
                    constraint.ub = constraint.lb + abs(r)
                elseif constraint_type == :L
                    constraint.lb = constraint.ub - abs(r)
                else
                    if r >= 0
                        constraint.ub = constraint.lb + r
                    else
                        constraint.lb = constraint.ub + r
                    end
                end
                ln = nextline(f)
            end
        end

        fill!(m.colLower, 0.0)

        if ln[1] == "BOUNDS"
            ln = nextline(f)
            while length(ln) > 1 
                
                boundType, varName = ln[1], ln[3]
                var = DicoVar[varName]
                if boundType == "LO"
                    val = parse(Float64, ln[4])
                    JuMP.setlowerbound(var, val)
                elseif boundType == "UP"
                    val = parse(Float64, ln[4])
                    JuMP.setupperbound(var, val)
                elseif boundType == "FX"
                    val = parse(Float64, ln[4])
                    JuMP.setupperbound(var, val)
                    JuMP.setlowerbound(var, val)
                elseif boundType == "FR"
                    JuMP.setupperbound(var, Inf)
                    JuMP.setlowerbound(var, -Inf)
                elseif boundType == "MI"
                    JuMP.setlowerbound(var, -Inf)
                elseif boundType == "PL"
                    JuMP.setupperbound(var, Inf)
                elseif boundType == "BV"
                    setcategory(var, :Bin)
                elseif boundType == "LI"
                    val = parse(Float64, ln[4])
                    setlowerbound(var, val)
                elseif boundType == "UI"
                    val = parse(Float64, ln[4])
                    setupperbound(var, val)
                elseif boundType == "SC"
                    error("Semi-continuous variables are not supported by vOptGeneric")
                end
                ln = nextline(f)
            end
        end
    end

    return m
end