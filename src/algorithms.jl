# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
function solve_lexico(m::JuMP.Model, verbose; args...)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E) #; empty!(vd.logObjs)
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
        solve_permutation(m, p, cstr_obj, cstr_rhs ; args...)
    end

    JuMP.delete.(m, cstr_obj)
    JuMP.delete.(m, cstr_rhs)
   
    return MOI.OPTIMAL
end

function solve_permutation(m::JuMP.Model, p, cstr_obj, cstr_rhs ; args...)

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

function solve_eps(m::JuMP.Model, ϵ::Float64, round_results, verbose ; args...)
    @info "ϵ = $ϵ"
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)# ; empty!(vd.logObjs)
    f1,f2 = vd.objs
    f1Sense,f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)
    
    #Set the first objective as an objective in the JuMP Model
    JuMP.set_objective(m, f1Sense, f1)
    
    R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
    R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
    weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

    #Declare the epsilon-constraint (RHS will be set later)
    RHS = JuMP.@variable(m)
    eps = (f2Sense == MOI.MIN_SENSE) ? JuMP.@constraint(m, f2 <= RHS) : JuMP.@constraint(m, f2 >= RHS)

    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    if status == MOI.OPTIMAL
        #While a solution exists
        while status == MOI.OPTIMAL
            #Get the score on the objectives
            f1Val = JuMP.value(f1)
            f2Val = JuMP.value(f2)
    
            #If last solution found is dominated by this one
            if length(vd.Y_N) > 0
                if weak_dom((f1Val, f2Val), vd.Y_N[end])
                    verbose && println(vd.Y_N[end], "dominated by ($f1Val, $f2Val)")
                    pop!(vd.Y_N) ; pop!(vd.X_E) #Remove last solution from Y_N and X_E
                end
            end

            #Store results in vOptData
            push!(vd.X_E, JuMP.value.(varArray))
            push!(vd.Y_N, [f1Val, f2Val])

            verbose && print("z1 = ", f1Val, ", z2 = ", f2Val)

            #Set the RHS of the epsilon-constraint
            if f2Sense == MOI.MIN_SENSE
                JuMP.fix(RHS, f2Val - ϵ)
                verbose && println(". Solving with f2 <= ", f2Val - ϵ)
            else
                JuMP.fix(RHS, f2Val + ϵ)
                verbose && println(". Solving with f2 >= ", f2Val + ϵ)
            end

            #And solve again
            JuMP.optimize!(m, ignore_optimize_hook=true)
            status = JuMP.termination_status(m)
        end

        #Sort X_E and Y_N
        s = sortperm(vd.Y_N, by = first)
        vd.X_E, vd.Y_N = vd.X_E[s], vd.Y_N[s]

        #Remove the epsilon constraint
        JuMP.delete(m, eps) ; JuMP.delete(m, RHS)
    else
        return status
    end
    return MOI.OPTIMAL
end

function solve_dicho(m::JuMP.Model, round_results, verbose; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)# ; empty!(vd.logObjs)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)

    #Set the first objective as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, f1Sense, f1)
    verbose && println("solving for z1")
    
    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    #If a solution exists
    if status == MOI.OPTIMAL

        yr_1 = JuMP.value(f1)
        yr_2 = JuMP.value(f2)

        #Store results in vOptData
        push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        push!(vd.X_E, JuMP.value.(varArray))
        # push!(vd.logObjs, f1)

        #Set the second objective as an objective in the JuMP JuMP.Model
        JuMP.set_objective(m, f2Sense, f2)
        verbose && println("solving for z2")

        #Solve with that objective
        JuMP.optimize!(m, ignore_optimize_hook=true)
        status = JuMP.termination_status(m)

        if status == MOI.OPTIMAL

            ys_1 = JuMP.value(f1)
            ys_2 = JuMP.value(f2)

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
                push!(vd.X_E, JuMP.value.(varArray))
                # push!(vd.logObjs, f2)
                dichoRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
            end
        
            #Sort X_E and Y_N
            s = sortperm(vd.Y_N, by = first)
            vd.Y_N = vd.Y_N[s]
            vd.X_E = vd.X_E[s]
            # vd.logObjs = vd.logObjs[s]

            R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
            R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
            weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

            #Filter X_E and Y_N :
            inds = Int[]
            for i = 1:length(vd.Y_N)-1
                if weak_dom(vd.Y_N[i], vd.Y_N[i+1])
                    push!(inds, i+1)
                elseif weak_dom(vd.Y_N[i+1], vd.Y_N[i])
                    push!(inds, i)
                end
            end
            deleteat!(vd.Y_N, inds)
            deleteat!(vd.X_E, inds)
            # deleteat!(vd.logObjs, inds)
        end
    end

    status
end

function dichoRecursion(m::JuMP.Model, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    f = AffExpr(0.0)

    if f1Sense==f2Sense
        lb = λ1*yr_1 + λ2*yr_2
        JuMP.set_objective(m, f1Sense, λ1*f1 + λ2*f2)
        verbose && println("solving for $λ1*f1 + $λ2*f2")    
        f = λ1*f1 + λ2*f2
    else
        lb = λ1*yr_1 - λ2*yr_2
        JuMP.set_objective(m, f1Sense, λ1*f1 - λ2*f2)
        verbose && println("solving for $λ1*f1 - $λ2*f2") 
        f = λ1*f1 - λ2*f2
    end

    JuMP.optimize!(m, ignore_optimize_hook=true)

    yt_1 = JuMP.value(f1)
    yt_2 = JuMP.value(f2)

    val = f1Sense == f2Sense ? λ1*yt_1 + λ2*yt_2 : λ1*yt_1 - λ2*yt_2

    if (f1Sense == MOI.MIN_SENSE && val < lb - 1e-4) || val > lb + 1e-4
        push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
        push!(vd.X_E, JuMP.value.(varArray))
        # push!(vd.logObjs, f)
        dichoRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, round_results, verbose ; args...)
        dichoRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
    end

end

function solve_Chalmet(m::JuMP.Model, step, verbose ; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)# ; empty!(vd.logObjs)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = JuMP.all_variables(m)

    #Set the first objective as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, f1Sense, f1)
    verbose && println("solving for z1")

    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    #If a solution exists
    if status == MOI.OPTIMAL

        yr_1 = JuMP.value(f1)
        yr_2 = JuMP.value(f2)

        #Store results in vOptData
        push!(vd.Y_N, [yr_1, yr_2])
        push!(vd.X_E, JuMP.value.(varArray))

        #Set the second objective as an objective in the JuMP JuMP.Model
        JuMP.set_objective(m, f2Sense, f2)
        verbose && println("solving for z2")


        #Solve with that objective
        JuMP.optimize!(m, ignore_optimize_hook=true)
        status = JuMP.termination_status(m)

        if status == MOI.OPTIMAL

            ys_1 = JuMP.value(f1)
            ys_2 = JuMP.value(f2)

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, [ys_1, ys_2])
                push!(vd.X_E, JuMP.value.(varArray))
                #Declare the constraints on z1 and z2 (RHS will be set later)
                rhs_z1 = JuMP.@variable(m)
                rhs_z2 = JuMP.@variable(m)
                cstr_z1 = f1Sense == MOI.MIN_SENSE ? JuMP.@constraint(m, f1 <= rhs_z1) : JuMP.@constraint(m, f1 >= rhs_z1)
                cstr_z2 = f2Sense == MOI.MIN_SENSE ? JuMP.@constraint(m, f2 <= rhs_z2) : JuMP.@constraint(m, f2 >= rhs_z2)

                ChalmetRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, rhs_z1, rhs_z2, step, verbose; args...)

                #Sort X_E and Y_N
                s = sortperm(vd.Y_N, by = first)
                vd.Y_N = vd.Y_N[s]
                vd.X_E = vd.X_E[s]

                R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
                R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
                weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

                #Filter X_E and Y_N :
                inds = Int[]
                for i = 1:length(vd.Y_N)-1
                    if weak_dom(vd.Y_N[i], vd.Y_N[i+1])
                        push!(inds, i+1)
                    elseif weak_dom(vd.Y_N[i+1], vd.Y_N[i])
                        push!(inds, i)
                    end
                end
                deleteat!(vd.Y_N, inds)
                deleteat!(vd.X_E, inds)

                JuMP.delete(m, cstr_z1)
                JuMP.delete(m, cstr_z2)
                JuMP.delete(m, rhs_z1)
                JuMP.delete(m, rhs_z2)
            end
        end
    end

    return MOI.OPTIMAL
end

function ChalmetRecursion(m::JuMP.Model, yr_1, yr_2, ys_1, ys_2, varArray, rhs_z1, rhs_z2, ϵ, verbose ; args...)

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    JuMP.set_objective(m, f1Sense, f1Sense==f2Sense ? f1 + f2 : f1 - f2)
    lbz1 = f1Sense==MOI.MAX_SENSE ? min(yr_1, ys_1) : max(yr_1, ys_1)
    lbz2 = f2Sense==MOI.MAX_SENSE ? min(yr_2, ys_2) : max(yr_2, ys_2)

    if f1Sense == MOI.MIN_SENSE
        JuMP.fix(rhs_z1, lbz1 - ϵ)
        verbose && print("solving with z1 <= ", lbz1 - ϵ)
    else
        JuMP.fix(rhs_z1, lbz1 + ϵ)
        verbose && print("solving with z1 >= ", lbz1 + ϵ)
    end

    if f2Sense == MOI.MIN_SENSE
        JuMP.fix(rhs_z2, lbz2 - ϵ)
        verbose && println(" and z2 <= ", lbz2 - ϵ)
    else
        JuMP.fix(rhs_z2, lbz2 + ϵ)
        verbose && println(" and z2 >= ", lbz2 + ϵ)
    end

    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    if status == MOI.OPTIMAL
        yt_1 = JuMP.value(f1)
        yt_2 = JuMP.value(f2)
        push!(vd.Y_N, [yt_1, yt_2])
        push!(vd.X_E, JuMP.value.(varArray))
        ChalmetRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, rhs_z1, rhs_z2, ϵ, verbose; args...)
        ChalmetRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, rhs_z1, rhs_z2, ϵ, verbose ; args...)
    end
end