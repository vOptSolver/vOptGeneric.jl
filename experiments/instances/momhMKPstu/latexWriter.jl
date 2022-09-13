using PyPlot
import PyPlot; const plt = PyPlot

function comparisonThreeMethods(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonTable.tex", "w")

    latex = raw"""\begin{table}[h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{Dichotomy}} & \multicolumn{2}{c}{\textbf{$\mathbf{\epsilon}$-constraint}}  & \multicolumn{3}{c}{\textbf{Branch-and-bound}} & \multicolumn{3}{c}{\textbf{Branch-and-cut}}
    \\
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-10} \cmidrule(r){11-13}
    ~ & ~ & ~ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$|\mathcal{X}_E|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$|\mathcal{X}_E|$} \\
    \midrule
    """
    println(fout, latex)
    methods = ["epsilon", "B&B", "B&C"] ; record_n = []
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods[2:end])

    for folder_n in readdir(work_dir * "/dicho") # ∀ filder_n
        count = 0
        avg_n = 0
        avg_m = 0
        avg_dicho_T = 0.0
        avg_dicho_Y = 0.0
        avg_ϵ_T = 0.0
        avg_ϵ_Y = 0.0
        avg_bb_T = 0.0
        avg_bb_Y = 0.0
        avg_bb_X = 0.0
        avg_bb_node = 0.0
        avg_bc_T = 0.0
        avg_bc_Y = 0.0
        avg_bc_X = 0.0
        avg_bc_node = 0.0

        for file in readdir(work_dir * "/dicho/" * string(folder_n) * "/") # ∀ file in dicho
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = []
            pts = [] ; X = []

            # write dichotomy result 
            include(work_dir * "/dicho/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")
            # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
            push!(times, total_times_used); push!(pts, size_Y_N)

            count += 1
            avg_n += vars
            avg_m += constr
            avg_dicho_T += total_times_used
            avg_dicho_Y += size_Y_N

            # write ϵ-constraint result (ϵ = 0.5 by default)
            if isfile(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
                include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N)

                avg_ϵ_T += total_times_used
                avg_ϵ_Y += size_Y_N
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
            end

            # write B&B result 
            if isfile(work_dir * "/bb/" * string(folder_n) * "/" * file)
                include(work_dir * "/bb/" * string(folder_n) * "/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N) 
                push!(X, size_X_E)

                avg_bb_T += total_times_used
                avg_bb_Y += size_Y_N
                avg_bb_X += size_X_E
                avg_bb_node += total_nodes
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
                push!(X, -1)
            end

            # write B&C result 
            if isfile(work_dir * "/bc/" * string(folder_n) * "/" * file)
                include(work_dir * "/bc/" * string(folder_n) * "/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N) 
                push!(X, size_X_E)

                avg_bc_T += total_times_used
                avg_bc_Y += size_Y_N
                avg_bc_X += size_X_E
                avg_bc_node += total_nodes
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
                push!(X, -1)
            end

            for i=1:4
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] == minimum(times)
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                elseif pts[i] == maximum(pts)
                    print(fout, " \\textbf{" * string(pts[i]) * "} & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end

                if i==3
                    X[1]==-1 ? print(fout, " - & ") : print(fout, string(X[1]) * " & ")
                elseif i==4
                    X[2]==-1 ? println(fout, " - \\\\") : println(fout, string(X[2]) * " \\\\")
                end
            end
    
        end

        avg_n = round(Int, avg_n/count)
        avg_m = round(Int, avg_m/count)
        avg_dicho_T = round(avg_dicho_T/count, digits = 2)
        avg_dicho_Y = round(avg_dicho_Y/count, digits = 2)
        avg_ϵ_T = round(avg_ϵ_T/count, digits = 2)
        avg_ϵ_Y = round(avg_ϵ_Y/count, digits = 2)
        avg_bb_T = round(avg_bb_T/count, digits = 2)
        avg_bb_Y = round(avg_bb_Y/count, digits = 2)
        avg_bb_X = round(avg_bb_X/count, digits = 2)
        avg_bb_node = round(avg_bb_node/count, digits = 2)
        avg_bc_T = round(avg_bc_T/count, digits = 2)
        avg_bc_Y = round(avg_bc_Y/count, digits = 2)
        avg_bc_X = round(avg_bc_X/count, digits = 2)
        avg_bc_node = round(avg_bc_node/count, digits = 2)

        append!(record_n, avg_n)
        append!(record_times["epsilon"], avg_ϵ_T) ; append!(record_times["B&B"], avg_bb_T) ; append!(record_times["B&C"], avg_bc_T) 
        append!(record_nodes["B&B"], avg_bb_node) ; append!(record_nodes["B&C"], avg_bc_node)


        println(fout, "\\cline{1-13} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" *
            string(avg_dicho_T) * "} & \\textbf{" * string(avg_dicho_Y) * "} & \\textbf{" * string(avg_ϵ_T) * "} & \\textbf{" *
            string(avg_ϵ_Y) * "} & \\textbf{" * string(avg_bb_T) * "} & \\textbf{" * string(avg_bb_Y) * "} & \\textbf{" * string(avg_bb_X) *
            "} & \\textbf{" * string(avg_bc_T) * "} & \\textbf{" * string(avg_bc_Y) * "} & \\textbf{" * string(avg_bc_X) * "} \\\\ \\cline{1-13}")

    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "\\label{tab:table_compare_$instances }")
    println(fout, "\\end{table}")
    close(fout)

    labels = [10, 20, 30, 40] ; loc = [0, 1, 2, 3] ; width = 0.3 # the width of the bars

    plt.bar(loc .- width, record_times["epsilon"], width, label="epsilon")
    plt.bar(loc, record_times["B&B"], width, label="B&B")
    plt.bar(loc .+ width, record_times["B&C"], width, label="B&C")
    plt.xticks(loc, labels)
    plt.xlabel("Number of variables")
    plt.ylabel("Computation time(s)", fontsize=14)
    plt.legend(methods)


    pos = loc .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["epsilon"][i]+0.5, s = record_times["epsilon"][i], size = 7)
    end

    pos = loc .+ (width) .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["B&B"][i]+0.5, s = record_times["B&B"][i], size = 7)
    end

    pos = loc .+ (2 * width) .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["B&C"][i]+0.5, s = record_times["B&C"][i], size = 7)
    end
    title("Influence of instance size on different algorithms' performance", fontsize=14)
    savefig(work_dir * "/comparisonTable.png")
    plt.close()

    # ---------------
    width = 0.4
    plt.bar(loc .- width/2, record_nodes["B&B"], width, label="B&B", color="darkorange")
    plt.bar(loc .+ width/2, record_nodes["B&C"], width, label="B&C", color="forestgreen")
    plt.xticks(loc, labels)
    plt.xlabel("Number of variables")
    plt.ylabel("Number of explored nodes", fontsize=14)
    plt.legend(methods[2:end])

    pos = loc .+ (width)/2# .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_nodes["B&B"][i]+0.5, s = record_nodes["B&B"][i], size = 7)
    end

    pos = loc .+ (3 *width)/2# .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_nodes["B&C"][i]+0.5, s = record_nodes["B&C"][i], size = 7)
    end
    title("Influence of instance size on different algorithms' tree size", fontsize=14)
    savefig(work_dir * "/comparisonNodes.png")
    plt.close()
end


function detailedMOBB_perform(instances::String)
    dir = "../../results/" * instances
    @assert isdir(dir) "This directory doesn't exist $dir !"

    fout = open(dir * "/detailedMOBB.tex", "w")

    latex = raw"""\begin{table}[!h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{4}{c}{\textbf{Time(s)}} & \multicolumn{2}{c}{\textbf{Nodes}}  & \textbf{Tree(MB)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$|\mathcal{X}_E|$}
    \\
    \cmidrule(r){4-7} \cmidrule(r){8-9} 
    ~ & ~ & ~ & \textbf{total} &\textbf{relax} & \textbf{dominance} & \textbf{incumbent} & \textbf{total} & \textbf{pruned} & ~ & ~ & ~ \\
    \midrule
    """
    println(fout, latex)
    labelT = ["BOLP", "dominance", "incumbent"] ; record_n = []
    labelNode = ["total", "pruned"]
    record_times = Dict(k => [] for k in labelT) ; record_nodes = Dict(k => [] for k in labelNode)

    for folder_n in readdir(dir * "/bb/")
        count = 0
        avg_n = 0
        avg_m = 0
        avg_totalT = 0.0
        avg_relaxT = 0.0
        avg_dominanceT = 0.0
        avg_incumbentT = 0.0
        avg_totalN= 0
        avg_prunedN = 0
        avg_treeSize = 0.0
        avg_YN = 0
        avg_XE = 0

        for file in readdir(dir * "/bb/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end
    
            print(fout, file * " & ")
            times = []
            pts = []
    
    
            include(dir * "/bb/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")
            print(fout, string(total_times_used)* " & " * string(relaxation_time) * " & " *
                string(test_dominance_time) * " & " * string(update_incumbent_time) * " & " *
                string(total_nodes) * " & " * string( pruned_nodes) * " & " * string(tree_size) * " & " *
                string(size_Y_N) * " & " * string(size_X_E)
            )
    
            println(fout, "\\\\")

            count += 1
            avg_n += vars
            avg_m += constr
            avg_totalT += total_times_used
            avg_relaxT += relaxation_time
            avg_dominanceT += test_dominance_time
            avg_incumbentT += update_incumbent_time
            avg_totalN += total_nodes
            avg_prunedN += pruned_nodes
            avg_treeSize += tree_size
            avg_YN += size_Y_N
            avg_XE += size_X_E
        end

        avg_n = round(Int, avg_n/count)
        avg_m = round(Int, avg_m/count)
        avg_totalT = round(avg_totalT/count, digits = 2)
        avg_relaxT = round(avg_relaxT/count, digits = 2)
        avg_dominanceT = round(avg_dominanceT/count, digits = 2)
        avg_incumbentT = round(avg_incumbentT/count, digits = 2)
        avg_totalN = round(avg_totalN/count, digits = 2)
        avg_prunedN = round(avg_prunedN/count, digits = 2)
        avg_treeSize = round(avg_treeSize/count, digits = 2)
        avg_YN = round(avg_YN/count, digits = 2)
        avg_XE = round(avg_XE/count, digits = 2)

        append!(record_n, avg_n)
        append!(record_times["BOLP"], avg_relaxT) ; append!(record_times["dominance"], avg_dominanceT) ; append!(record_times["incumbent"], avg_incumbentT)
        append!(record_nodes["total"], avg_totalN) ; append!(record_nodes["pruned"], avg_prunedN)

        println(fout, "\\cline{1-12} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" *
            string(avg_totalT) * "} & \\textbf{" * string(avg_relaxT) * "} & \\textbf{" * string(avg_dominanceT) * "} & \\textbf{" *
            string(avg_incumbentT) * "} & \\textbf{" * string(avg_totalN) * "} & \\textbf{" * string(avg_prunedN) * "} & \\textbf{" *
            string(avg_treeSize) * "} & \\textbf{" * string(avg_YN) * "} & \\textbf{" * string(avg_XE) *"} \\\\ \\cline{1-12}")
    end

    latex = raw"""\bottomrule
    \end{tabular}%
    }%
    \caption{The detailed experimental information about BO01B\&B algorithm.}
    \label{tab:table_bb}
    \end{table}
    """
    println(fout, latex)
    close(fout)

    labels = [10, 20, 30, 40] ; loc = [0, 1, 2, 3] ; width = 0.3 # the width of the bars

    plt.bar(loc .- width, record_times["BOLP"], width, label="BOLP")
    plt.bar(loc, record_times["dominance"], width, label="dominance")
    plt.bar(loc .+ width, record_times["incumbent"], width, label="incumbent")
    plt.xticks(loc, labels)
    plt.xlabel("Number of variables", fontsize=14)
    plt.ylabel("Computation time(s)", fontsize=14)
    plt.legend(labelT)

    pos = loc .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["BOLP"][i]+0.5, s = record_times["BOLP"][i], size = 7)
    end

    pos = loc .+ (width) .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["dominance"][i]+0.5, s = record_times["dominance"][i], size = 7)
    end

    pos = loc .+ (2 * width) .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["incumbent"][i]+0.5, s = record_times["incumbent"][i], size = 7)
    end
    title("The influence of instance size on BO01B&B algorithm", fontsize=14)
    savefig(dir * "/BBperformTimes.png")
    plt.close()


    width = 0.4
    plt.bar(loc .- width/2, record_nodes["total"], width, label="total")
    plt.bar(loc .+ width/2, record_nodes["pruned"], width, label="pruned")
    plt.xticks(loc, labels)
    plt.xlabel("Number of variables", fontsize=14)
    plt.ylabel("Number of nodes", fontsize=14)
    plt.legend(labelNode)

    pos = loc .+ (width)/2# .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_nodes["total"][i]+0.5, s = record_nodes["total"][i], size = 7)
    end

    pos = loc .+ (3 *width)/2# .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_nodes["pruned"][i]+0.5, s = record_nodes["pruned"][i], size = 7)
    end
    title("The influence of instance size on B&B tree", fontsize=14)
    savefig(dir * "/BBperformNodes.png")
    plt.close()
end

function MOBC_perform(instances::String)
    bc = "/bc"
    dir = "../../results/" * instances * bc
    @assert isdir(dir) "This directory doesn't exist $dir !"

    fout = open(dir * "/MOBC_bc.tex", "w")

    latex = raw"""\begin{sidewaystable}[h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{Nodes}} & \multicolumn{2}{c}{\textbf{CP iterations}} & \multicolumn{3}{c}{\textbf{Cuts applied}} & \textbf{Cuts} & \multicolumn{5}{c}{\textbf{CP Time(s)}} & \textbf{B\&C Time(s)} & \textbf{$|\mathcal{Y}_N|$}
    \\
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-10} \cmidrule(r){12-16}
    ~ & ~ & ~ & \textbf{total} & \textbf{pruned} & \textbf{total} & \textbf{average} & \textbf{total} & \textbf{sp} & \textbf{mp} & ~ & \textbf{total} & \textbf{dichotomy} & \textbf{pool oper} & \textbf{separators} & \textbf{cuts oper} & ~ & ~ \\
    \midrule
    """
    println(fout, latex)

    for folder_n in readdir(dir)
        if !isdir(dir * "/" * string(folder_n) ) continue end 
        count = 0
        avg_n = 0 ; avg_m = 0
        avg_totalNodes= 0 ; avg_prunedNodes = 0
        avg_totalIte = 0 ; avg_pernodesIte = 0
        avg_totalCuts = 0 ; avg_spCuts = 0 ; avg_mpCuts = 0 ; avg_Cuts = 0
        avg_totalT = 0.0 ; avg_dichoT = 0.0 ; avg_poolT = 0.0 ; avg_sepaT = 0.0 ; avg_cutsT = 0.0
        avg_BCtime = 0.0 ; avg_YN = 0

        for file in readdir(dir * "/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end
    
            print(fout, file * " & ")
    
            include(dir * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")
            print(fout, string(total_nodes) * " & " * string(pruned_nodes) * " & " * string(ite_total) * " & "*
                string( round(ite_total/total_nodes, digits = 2) ) * " & " * string(cuts_applied) * " & " * string(sp_cuts)*
                " & " * string(mp_cuts) * " & " * string(cuts_total) * " & " * string(times_total_for_cuts) * " & " *
                string(times_calling_dicho) * " & " * string(times_oper_cutPool) * " & " *string(times_calling_separators) * " & "*
                string(times_add_retrieve_cuts) * " & " *string(total_times_used) * " & " *string(size_Y_N)
            )

            println(fout, "\\\\")

            count += 1
            avg_n += vars ; avg_m += constr
            avg_totalNodes += total_nodes ; avg_prunedNodes += pruned_nodes
            avg_totalIte += ite_total ; avg_pernodesIte += round(ite_total/total_nodes, digits = 2)
            avg_totalCuts += cuts_applied ; avg_spCuts += sp_cuts ; avg_mpCuts += mp_cuts ; avg_Cuts += cuts_total
            avg_totalT += times_total_for_cuts ; avg_dichoT += times_calling_dicho ; avg_poolT += times_oper_cutPool ; avg_sepaT += times_calling_separators ; avg_cutsT += times_add_retrieve_cuts
            avg_BCtime += total_times_used ; avg_YN += size_Y_N
        end

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        avg_totalNodes = round(avg_totalNodes/count, digits = 2)
        avg_prunedNodes = round(avg_prunedNodes/count, digits = 2)
        avg_totalIte = round(avg_totalIte/count, digits = 2)
        avg_pernodesIte = round(avg_pernodesIte/count, digits = 2)
        avg_totalCuts = round(avg_totalCuts/count, digits = 2)
        avg_spCuts = round(avg_spCuts/count, digits = 2)
        avg_mpCuts = round(avg_mpCuts/count, digits = 2)
        avg_Cuts = round(avg_Cuts/count, digits = 2)
        avg_dichoT = round(avg_dichoT/count, digits = 2)
        avg_totalT = round(avg_totalT/count, digits = 2)
        avg_poolT = round(avg_poolT/count, digits = 2)
        avg_sepaT = round(avg_sepaT/count, digits = 2)
        avg_cutsT = round(avg_cutsT/count, digits = 2)
        avg_BCtime = round(avg_BCtime/count, digits = 2)
        avg_YN = round(avg_YN/count, digits = 2)


        println(fout, "\\cline{1-18} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" *
            string(avg_totalNodes) * "} & \\textbf{" * string(avg_prunedNodes) * "} & \\textbf{" * string(avg_totalIte) * "} & \\textbf{" *
            string(avg_pernodesIte) * "} & \\textbf{" * string(avg_totalCuts) * "} & \\textbf{" * string(avg_spCuts) * "} & \\textbf{" *
            string(avg_mpCuts) * "} & \\textbf{" * string(avg_Cuts) * "} & \\textbf{" * string(avg_dichoT) *"} & \\textbf{" *
            string(avg_totalT) * "} & \\textbf{" * string(avg_poolT) * "} & \\textbf{" * string(avg_sepaT) * "} & \\textbf{" *
            string(avg_cutsT) * "} & \\textbf{" * string(avg_BCtime) * "} & \\textbf{" * string(avg_YN) * "} " * "\\\\ \\cline{1-18}")
    end


    latex = raw"""\bottomrule
    \end{tabular}%
    }%
    \caption{.}
    \label{tab:table_bc}
    \end{sidewaystable}
    """
    println(fout, latex)
    close(fout)

end



detailedMOBB_perform("momhMKPstu/MOBKP/set3")

comparisonThreeMethods("momhMKPstu/MOBKP/set3")

MOBC_perform("momhMKPstu/MOBKP/set3")