# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.

function solve_eps(m::JuMP.Model; step::Float64=1., round_results::Bool=false, verbose, kwargs...)
    ϵ = step
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
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
