function solve_eps(m::Model, ϵ::Float64)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    #Set the first objective as an objective in the JuMP model
    m.obj=f1
    m.objSense = f1Sense

    #Solve with that objective
    status = solve(m, ignore_solve_hook=true)

    #If a solution exists
    if status == :Optimal

        #Declare the epsilon-constraint (RHS will be set later)
        if f2Sense == :Min        
            eps = @constraint(m, f2.aff <= 0.0)
        else
            eps = @constraint(m, f2.aff >= 0.0)
        end

        varArray = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

        #While a solution exists
        while status == :Optimal
            #Get the score on the objectives
            f1Val = m.objVal
            f2Val = JuMP.getvalue(f2)

            R1 = f1Sense==:Min ? (<=) : (>=)
            R2 = f2Sense==:Min ? (<=) : (>=)

            #If last solution found is dominated by this one
            if !isempty(vd.Y_N)
                Y_m1 = vd.Y_N[end]
                if R1(f1Val,Y_m1[1]) && R2(f2Val, Y_m1[2])
                    #Remove last solution from Y_N and X_E
                    pop!(vd.Y_N) ; pop!(vd.X_E)
                end
            end

            #If this solution isn't dominated by the last one :
            if isempty(vd.Y_N) || !(R1(vd.Y_N[end][1],f1Val) && R2(vd.Y_N[end][1], f2Val))
                #Store results in vOptData
                push!(vd.Y_N, (f1Val, f2Val))
                push!(vd.X_E, JuMP.getvalue.(varArray))
            end

            print("z1 = ", f1Val, ", z2 = ", f2Val)
            #Set the RHS of the epsilon-constraint
            if f2Sense == :Min
                JuMP.setRHS(eps, f2Val - f2.aff.constant - ϵ)
                println(". Solving with f2 <= ", f2Val - ϵ)
            else
                JuMP.setRHS(eps, f2Val - f2.aff.constant + ϵ)
                println(". Solving with f2 >= ", f2Val + ϵ)
            end

            # for (k,v) in varDict
            #     setvalue(v, getvalue(v))
            # end

            #And solve again
            status = solve(m, ignore_solve_hook=true, suppress_warnings=true)
        end
        ###
        #JuMP doesn't support removing constraints
        #To leave the model unaltered we change the value of the RHS
        ###
        if f2Sense == :Min
            JuMP.setRHS(eps, Float64(typemax(Int)))
        else
            JuMP.setRHS(eps, Float64(typemin(Int)))
        end
    else
        return status
    end
    return :Optimal
end

function solve_dicho(m::Model)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    #Set the first objective as an objective in the JuMP model
    m.obj=f1
    m.objSense = f1Sense

    #Solve with that objective
    status = solve(m, ignore_solve_hook=true)

    #If a solution exists
    if status == :Optimal

        yr_1 = m.objVal
        yr_2 = JuMP.getvalue(f2)

        #Store results in vOptData
        push!(vd.Y_N, (yr_1, yr_2))
        push!(vd.X_E, JuMP.getvalue.(varArray))


        #Set the second objective as an objective in the JuMP model
        m.obj=f2
        m.objSense = f2Sense

        #Solve with that objective
        status = solve(m, ignore_solve_hook=true)

        if status == :Optimal

            ys_1 = JuMP.getvalue(f1)
            ys_2 = m.objVal

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, (ys_1, ys_2))
                push!(vd.X_E, JuMP.getvalue.(varArray))
                dichoRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray)
            end
        
            #Lazy sort X_E and Y_N
            s = sortperm(vd.Y_N, by = e -> e[1])
            vd.Y_N = vd.Y_N[s]
            vd.X_E = vd.X_E[s]

            #We can use weak dominance definition here
            dom_min_min(a,b) = a[1] <= b[1] && a[2] <= b[2]
            dom_max_max(a,b) = a[1] >= b[1] && a[2] >= b[2]
            dom_min_max(a,b) = a[1] <= b[1] && a[2] >= b[2]
            dom_max_min(a,b) = a[1] >= b[1] && a[2] <= b[2]

            if f1Sense==f2Sense==:Min
                dominates = dom_min_min
            elseif f1Sense==f2Sense==:Max
                dominates = dom_max_max
            elseif f1Sense==:Min
                dominates = dom_min_max
            else
                dominates = dom_max_min
            end

            #Filter X_E and Y_N :
            inds = Int[]
            for i = 1:length(vd.Y_N)-1
                if dominates(vd.Y_N[i], vd.Y_N[i+1])
                    push!(inds, i+1)
                elseif dominates(vd.Y_N[i+1], vd.Y_N[i])
                    push!(inds, i)
                end
            end
            deleteat!(vd.Y_N, inds)
            deleteat!(vd.X_E, inds)
        end
    end

    status
end

function dichoRecursion(m::Model, yr_1, yr_2, ys_1, ys_2, varArray)

    vd = getvOptData(m)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    m.objSense = f1Sense

    if f1Sense==f2Sense
        lb = λ1*yr_1 + λ2*yr_2
        m.obj = λ1*f1 + λ2*f2
    else
        lb = λ1*yr_1 - λ2*yr_2
        m.obj = λ1*f1 - λ2*f2
    end

    solve(m, ignore_solve_hook=true)

    yt_1 = JuMP.getvalue(f1)
    yt_2 = JuMP.getvalue(f2)


    if f1Sense==f2Sense
        val = λ1*yt_1 + λ2*yt_2
    else
        val = λ1*yt_1 - λ2*yt_2
    end

    if f1Sense == :Min && val < lb - 1e-4 || val > lb + 1e-4
        push!(vd.Y_N, (yt_1, yt_2))
        push!(vd.X_E, JuMP.getvalue.(varArray))
        dichoRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray)
        dichoRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray)
    end

end


function solve_Chalmet(m::Model, step)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    #Set the first objective as an objective in the JuMP model
    m.obj=f1
    m.objSense = f1Sense

    #Solve with that objective
    status = solve(m, ignore_solve_hook=true)

    #If a solution exists
    if status == :Optimal

        yr_1 = m.objVal
        yr_2 = JuMP.getvalue(f2)

        #Store results in vOptData
        push!(vd.Y_N, (yr_1, yr_2))
        push!(vd.X_E, JuMP.getvalue.(varArray))


        #Set the second objective as an objective in the JuMP model
        m.obj=f2
        m.objSense = f2Sense

        #Solve with that objective
        status = solve(m, ignore_solve_hook=true)

        if status == :Optimal

            ys_1 = JuMP.getvalue(f1)
            ys_2 = m.objVal

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, (ys_1, ys_2))
                push!(vd.X_E, JuMP.getvalue.(varArray))
                #Declare the constraints on z1 and z2 (RHS will be set later)
                if f1Sense == :Min
                    cstrz1 = @constraint(m, f1.aff <= 0.0)
                else
                    cstrz1 = @constraint(m, f1.aff >= 0.0)
                end
                if f2Sense == :Min
                    cstrz2 = @constraint(m, f2.aff <= 0.0)
                else
                    cstrz2 = @constraint(m, f2.aff >= 0.0)
                end
                ChalmetRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, cstrz1, cstrz2, step)

                #Lazy sort X_E and Y_N
                s = sortperm(vd.Y_N, by = e -> e[1])
                vd.Y_N = vd.Y_N[s]
                vd.X_E = vd.X_E[s]

                #We can use weak dominance definition here
                dom_min_min(a,b) = a[1] <= b[1] && a[2] <= b[2]
                dom_max_max(a,b) = a[1] >= b[1] && a[2] >= b[2]
                dom_min_max(a,b) = a[1] <= b[1] && a[2] >= b[2]
                dom_max_min(a,b) = a[1] >= b[1] && a[2] <= b[2]

                if f1Sense==f2Sense==:Min
                    dominates = dom_min_min
                elseif f1Sense==f2Sense==:Max
                    dominates = dom_max_max
                elseif f1Sense==:Min
                    dominates = dom_min_max
                else
                    dominates = dom_max_min
                end

                #Filter X_E and Y_N :
                inds = Int[]
                for i = 1:length(vd.Y_N)-1
                    if dominates(vd.Y_N[i], vd.Y_N[i+1])
                        push!(inds, i+1)
                    elseif dominates(vd.Y_N[i+1], vd.Y_N[i])
                        push!(inds, i)
                    end
                end
                deleteat!(vd.Y_N, inds)
                deleteat!(vd.X_E, inds)
            end
        end
    end

    return status
end

function ChalmetRecursion(m::Model, yr_1, yr_2, ys_1, ys_2, varArray, cstrz1, cstrz2, ϵ)

    vd = getvOptData(m)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    m.objSense = f1Sense
    if f1Sense==f2Sense
        m.obj = f1 + f2
    else
        m.obj = f1 - f2
    end

    lbz1 = f1Sense==:Max ? min(yr_1, ys_1) : max(yr_1, ys_1)
    lbz2 = f2Sense==:Max ? min(yr_2, ys_2) : max(yr_2, ys_2)

    if f1Sense == :Min
        JuMP.setRHS(cstrz1, lbz1 - f1.aff.constant - ϵ)
        print(". Solving with f1 <= ", lbz1 - ϵ)
    else
        JuMP.setRHS(cstrz1, lbz1 - f1.aff.constant + ϵ)
        print(". Solving with f1 >= ", lbz1 + ϵ)
    end

    if f2Sense == :Min
        JuMP.setRHS(cstrz2, lbz2 - f2.aff.constant - ϵ)
        println(" and f2 <= ", lbz2 - ϵ)
    else
        JuMP.setRHS(cstrz2, lbz2 - f2.aff.constant + ϵ)
        println(" and f2 >= ", lbz2 + ϵ)
    end

    status = solve(m, ignore_solve_hook=true)

    if status == :Optimal

        yt_1 = JuMP.getvalue(f1)
        yt_2 = JuMP.getvalue(f2)

        push!(vd.Y_N, (yt_1, yt_2))
        push!(vd.X_E, JuMP.getvalue.(varArray))
        ChalmetRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, cstrz1, cstrz2, ϵ)
        ChalmetRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, cstrz1, cstrz2, ϵ)

    end
end