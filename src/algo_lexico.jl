# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.

function solve_lexico(m::JuMP.Model, verbose; kwargs...)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    objs = vd.objs
    objSenses = vd.objSenses
    nbObj = length(objs)

    #Check that problem is feasible and that no objective is unbounded
    for i = 1:nbObj
        JuMP.set_objective(m, objSenses[i], objs[i])
        JuMP.optimize!(m, ignore_optimize_hook=true)
        status = JuMP.termination_status(m)
        if status != MOI.OPTIMAL
            return status
        end
    end

    #Create the constraints on the objectives
    cstr_rhs = [JuMP.@variable(m) for _ = 1:nbObj]
    cstr_obj = [objSenses[i] == MOI.MAX_SENSE ? JuMP.@constraint(m, objs[i] >= cstr_rhs[i]) : JuMP.@constraint(m, objs[i] <= cstr_rhs[i]) for i=1:nbObj]

    for p in Combinatorics.permutations(1:nbObj, nbObj)
        verbose && println("solving for objectives $p")
        solve_permutation(m, p, cstr_obj, cstr_rhs ; kwargs...)
    end

    JuMP.delete.(m, cstr_obj)
    JuMP.delete.(m, cstr_rhs)

    return MOI.OPTIMAL
end

function solve_permutation(m::JuMP.Model, p, cstr_obj, cstr_rhs ; kwargs...)

    vd = getvOptData(m)
    objs = vd.objs
    objSenses = vd.objSenses

    #Set the first objective of the permutation as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, objSenses[first(p)], objs[first(p)])

    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)

    status = JuMP.termination_status(m)
    if status != MOI.OPTIMAL
        return status
    end

    for i = 2:length(p)
        fVal = JuMP.value(objs[p[i-1]]) #get the value for the last objective solved
        slack = objSenses[p[i-1]] == MOI.MAX_SENSE ? -1e-8 : 1e-8
        # JuMP.fix(cstr_rhs[p[i-1]], fVal - objs[p[i-1]].constant + slack)
        JuMP.fix(cstr_rhs[p[i-1]], fVal + slack)
        JuMP.set_objective(m, objSenses[p[i]], objs[p[i]]) #set the i-th objective of the permutation in the JuMP JuMP.Model
        JuMP.optimize!(m, ignore_optimize_hook=true) #and solve
    end

    varArray = JuMP.all_variables(m)
    #Store results in vOptData
    push!(vd.Y_N, map(JuMP.value, objs))
    push!(vd.X_E, JuMP.value.(varArray))

    for var in cstr_rhs
        JuMP.is_fixed(var) && JuMP.unfix.(var) #errors if calling unfix on an unfixed variable
    end

end
