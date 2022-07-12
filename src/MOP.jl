# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.

function _varname(v::JuMP.VariableRef)
    name = JuMP.name(v)
    for (pat, sub) in [("[", "_"), ("]", ""), (",", "_")]
        name = replace(name, pat => sub)
    end
    name
end

function writeMOP(m, fname::AbstractString, genericnames = false)

    open(fname, "w") do f

        print_short = x -> Base.Grisu.print_shortest(f, x)

        print(f, "NAME   vOptModel\n")
    
        cstrEqualTo, cstrGreaterThan, cstrLessThan, cstrInterval = [JuMP.all_constraints(m, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, set_type) 
            for set_type in (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64})]
        numRows = sum(length, (cstrEqualTo, cstrGreaterThan, cstrLessThan, cstrInterval))

        md = getvOptData(m)

        print(f, "OBJSENSE\n")
        if md.objSenses[1] == MOI.MIN_SENSE
            print(f, " MIN\n")
        else
            print(f, " MAX\n")
        end

        flippedObj = [sense != md.objSenses[1] for sense in md.objSenses]
        if any(flippedObj)
            @info "MOP format does not support different optimisation senses. Flipping objective coefficients."
        end

        # Objective and constraint names
        print(f, "ROWS\n")
        for i = 1:length(md.objs)
            print(f, " N  OBJ$i\n")
        end

        dictCstrIndex = Dict{JuMP.ConstraintRef, Int}()
        cstr_index = 1
        for cstrRef in cstrLessThan
            print(f, " L  CON$cstr_index\n")
            push!(dictCstrIndex, cstrRef => cstr_index)
            cstr_index += 1
        end
        for cstrRef in cstrEqualTo
            print(f, " E  CON$cstr_index\n")
            push!(dictCstrIndex, cstrRef => cstr_index)
            cstr_index += 1
        end
        for cstrRef in cstrGreaterThan
            print(f, " G  CON$cstr_index\n")
            push!(dictCstrIndex, cstrRef => cstr_index)
            cstr_index += 1
        end
        for cstrRef in cstrInterval
            print(f, " E  CON$cstr_index\n")
            push!(dictCstrIndex, cstrRef => cstr_index)
            cstr_index += 1
        end

        # Output each column
        inintegergroup = false
        print(f, "COLUMNS\n")
        for var in JuMP.all_variables(m)
            if JuMP.is_binary(var) || JuMP.is_integer(var) && !inintegergroup
                print(f, "    MARKER    'MARKER'                 'INTORG'\n")
                inintegergroup = true
            elseif !JuMP.is_binary(var) && !JuMP.is_integer(var) && inintegergroup
                print(f, "    MARKER    'MARKER'                 'INTEND'\n")
                inintegergroup = false
            end

            for set_type in (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64})
                for cstrRef in JuMP.all_constraints(m, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, set_type)
                    con = JuMP.constraint_object(cstrRef)
                    terms = JuMP.linear_terms(con.func)
                    for (coeff, term) in terms
                        if var == term
                            print(f, "    $(_varname(var))  CON$(dictCstrIndex[cstrRef])  ") ; print_short(coeff) ; println(f)
                        end
                    end
                end
            end

            for i = 1:length(md.objs)
                for (coeff, term) in JuMP.linear_terms(md.objs[i])
                    if var == term
                        print(f,"    $(_varname(var))  OBJ$i  ") ; print_short(flippedObj[i] ? -coeff : coeff) ; println(f)
                    end
                end
            end

        end
        if inintegergroup
            print(f, "    MARKER    'MARKER'                 'INTEND'\n")
        end

        # RHSs
        print(f, "RHS\n")

        # TODO : write objectives constant

        cstr_index = 1
        for cstrRef in cstrLessThan
            con = JuMP.constraint_object(cstrRef)
            rhs = con.set.upper
            print(f, "    RHS    CON$cstr_index    ") ; print_short(rhs) ; println(f)
            cstr_index += 1
        end
        for cstrRef in cstrEqualTo
            con = JuMP.constraint_object(cstrRef)
            rhs = con.set.value
            print(f, "    RHS    CON$cstr_index    ") ; print_short(rhs) ; println(f)
            cstr_index += 1
        end
        for cstrRef in cstrGreaterThan
            con = JuMP.constraint_object(cstrRef)
            rhs = con.set.lower
            print(f, "    RHS    CON$cstr_index    ") ; print_short(rhs) ; println(f)
            cstr_index += 1
        end
        for cstrRef in cstrInterval
            con = JuMP.constraint_object(cstrRef)
            rhs = con.set.lower
            print(f, "    RHS    CON$cstr_index    ") ; print_short(rhs) ; println(f)
            cstr_index += 1
        end

        # RANGES
        if length(cstrInterval) > 0
            print(f, "RANGES\n")
            cstr_index = numRows - length(cstrInterval)
            for cstrRef in cstrInterval
                con = JuMP.constraint_object(cstrRef)
                lb = con.set.lower
                ub = con.set.upper
                print(f, "    rhs    CON$cstr_index    ") ; print_short(ub - lb) ; println(f)
            end
        end

        # BOUNDS
        print(f, "BOUNDS\n")
        for var in JuMP.all_variables(m)
            if JuMP.is_binary(var)
                #Binary
                print(f, "  BV BOUND $(_varname(var))") ; println(f)
            elseif JuMP.is_fixed(var)
                # Fixed variable
                print(f, "  FX BOUND $(_varname(var)) ") ; print_short(JuMP.fix_value(var)) ; println(f)
            elseif JuMP.has_lower_bound(var) && JuMP.lower_bound(var) == 0 && JuMP.has_upper_bound(var) && JuMP.upper_bound(var) != Inf
                # Default lower 0, and an upper
                print(f, "  UP BOUND $(_varname(var))") ; print_short(JuMP.upper_bound(var)) ; println(f)
            elseif JuMP.has_lower_bound(var) && JuMP.lower_bound(var) == 0 && (!JuMP.has_upper_bound(var) || JuMP.upper_bound(var) == Inf)
                # Default bounds.
                print(f, "  PL BOUND $(_varname(var))")
                println(f)
            elseif (!JuMP.has_lower_bound(var) || JuMP.lower_bound(var) == -Inf) && (!JuMP.has_upper_bound(var) || JuMP.upper_bound(var) == Inf)
                # Free
                print(f, "  FR BOUND $(_varname(var))\n")
            elseif JuMP.has_lower_bound(var) && JuMP.lower_bound(var) != -Inf && (!JuMP.has_upper_bound(var) || JuMP.upper_bound(var) == +Inf)
                # No upper, but a lower
                print(f, "  PL BOUND $(_varname(var))\n")
                print(f, "  LO BOUND $(_varname(var))") ; print_short(JuMP.lower_bound(var)) ; println(f)
            elseif (!JuMP.has_lower_bound(var) || JuMP.lower_bound(var) == -Inf) && JuMP.has_upper_bound(var) && JuMP.upper_bound(var) != +Inf
                # No lower, but a upper
                print(f, "  MI BOUND $(_varname(var))\n")
                print(f, "  UP BOUND $(_varname(var)) ") ; print_short(JuMP.upper_bound(var)) ; println(f)
            else
                # Lower and upper
                print(f, "  LO BOUND $(_varname(var)) ") ; print_short(JuMP.lower_bound(var)) ; println(f)
                print(f, "  UP BOUND $(_varname(var)) ") ; print_short(JuMP.upper_bound(var)) ; println(f)
            end
        end
    
        print(f, "ENDATA\n")
    end
    nothing
end

function parseMOP(fname::AbstractString, optimizer_factory = nothing)
    
    m = vModel(optimizer_factory) ; JuMP.set_silent(m)
    nextline = (f) -> split(chomp(readline(f)), ' ', keepempty = false)

    open(fname) do f
        ln = nextline(f)
        while ln[1] != "NAME"
            ln = nextline(f)
        end

        objSense = MOI.MIN_SENSE
        ln = nextline(f)
        if ln[1] == "OBJSENSE"
            ln = nextline(f)
            if ln[1] == "MAXIMIZE" || ln[1] == "MAX"
                objSense = MOI.MAX_SENSE
            end
            ln = nextline(f)
        end

        if ln[1] == "ROWS"
            ln = nextline(f)
            dictExpr = Dict{String,Tuple{JuMP.GenericAffExpr,Symbol}}()
            objOrder = Vector{String}()
            while ln[1] != "COLUMNS"
                push!(dictExpr, ln[2] => (JuMP.AffExpr(), Symbol(ln[1])))
                ln[1] == "N" && push!(objOrder, ln[2])
                ln = nextline(f)
            end
        end

        if ln[1] == "COLUMNS"
            ln = nextline(f)
            dictVar = Dict{String,JuMP.VariableRef}()
            isInt = false
            while ln[1] != "RHS"
                if ln[3] == "'INTORG'"
                    isInt = true
                elseif ln[3] == "'INTEND'"
                    isInt = false
                else
                    if !(ln[1] in keys(dictVar))
                        v = JuMP.@variable(m, integer = isInt, base_name = String(ln[1]), lower_bound = 0.0)
                        push!(dictVar, ln[1] => v)
                    end
                    expr, expr_type = dictExpr[ln[2]]
                    var = dictVar[ln[1]]
                    coeff = parse(Float64, ln[3])

                    dictExpr[ln[2]] = (expr + coeff * var, expr_type)
                end
                ln = nextline(f)
            end
        end

        dictObj = Dict{String,JuMP.GenericAffExpr}()
        for (k, (expr, expr_type)) in dictExpr
            if expr_type == :N
                push!(dictObj, k => expr)
            end
        end
        filter!(p->p |> last |> last != :N, dictExpr)
        
        md = m.ext[:vOpt]
        for key in objOrder
            push!(md.objs, dictObj[key])
            push!(md.objSenses, objSense)
        end
        
        # TODO : read objectives constant
        dictRHS = Dict{String,Float64}()
        if ln[1] == "RHS"
            ln = nextline(f)
            while length(ln) > 1
                cstr_name = ln[2]
                rhs = parse(Float64, ln[3])
                # cstr_name == "nonS@800" && error()
                push!(dictRHS, cstr_name => rhs)
                ln = nextline(f)
            end
        end

        if ln[1] == "RANGES"
            ln = nextline(f)
            while length(ln) > 1
                cstr_name = ln[2]
                expr, constraint_type = dictExpr[cstr_name]
                rhs = dictRHS[cstr_name]
                r = parse(Float64, ln[3])
                if constraint_type == :G
                    JuMP.@constraint(m, rhs <= expr <= rhs + abs(r))
                elseif constraint_type == :L
                    JuMP.@constraint(m, rhs - abs(r) <= expr <= rhs)
                else
                    if r >= 0
                        JuMP.@constraint(m, rhs <= expr <= rhs + r)                        
                    else
                        JuMP.@constraint(m, rhs + r <= expr <= rhs)                        
                    end
                end
                delete!(dictExpr, cstr_name)
                delete!(dictRHS, cstr_name)
                ln = nextline(f)
            end
        end

        # add all constraints that aren't range constraints
        for (cstr_name, rhs) in dictRHS
            expr, constraint_type = dictExpr[cstr_name]
            if constraint_type == :E
                JuMP.@constraint(m, expr == rhs)
            elseif constraint_type == :G
                JuMP.@constraint(m, expr >= rhs)
            elseif constraint_type == :L
                JuMP.@constraint(m, expr <= rhs)
            end
            delete!(dictExpr, cstr_name)
        end

        undeclared_variables = String[]
        if ln[1] == "BOUNDS"
            ln = nextline(f)
            while length(ln) > 1 
                boundType, varName = ln[1], ln[3]
                val = length(ln) > 3 ? parse(Float64, ln[4]) : 0.
                if haskey(dictVar, varName) 
                    var = dictVar[varName]
                    if boundType == "LO"
                        JuMP.set_lower_bound(var, val)
                    elseif boundType == "UP"
                        JuMP.set_upper_bound(var, val)
                    elseif boundType == "FX"
                        JuMP.fix(var, val)
                    elseif boundType == "FR"
                        JuMP.set_upper_bound(var, Inf)
                        JuMP.set_lower_bound(var, -Inf)
                    elseif boundType == "MI"
                        JuMP.set_lower_bound(var, -Inf)
                    elseif boundType == "PL"
                        JuMP.set_upper_bound(var, Inf)
                    elseif boundType == "BV"
                        JuMP.is_integer(var) && JuMP.unset_integer(var)
                        JuMP.set_binary(var)
                    elseif boundType == "LI"
                        JuMP.set_lower_bound(var, val)
                    elseif boundType == "UI"
                        JuMP.set_upper_bound(var, val)
                    elseif boundType == "SC"
                        error("Semi-continuous variables are not supported by vOptGeneric")
                    else
                        error("Unrecognised bound type : $boundType")
                    end
                else 
                    push!(undeclared_variables, varName)
                end
                ln = nextline(f)
            end
        end
        isempty(undeclared_variables) || @warn "Some variables were declared in BOUNDS but were never used, they were not added to the model" variables = undeclared_variables 
        isempty(dictExpr) || @warn "Some constraints are missing a RHS and were not added to the model" constraints = keys(dictExpr) |> collect
    end

    return m
end