using PyPlot
import PyPlot; const plt = PyPlot

function comparisonMethods(instances::String)
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

    for file in readdir(work_dir * "/dicho/") # ∀ file in dicho
        if split(file, ".")[end] == "png"
            continue
        end

        instName = split(file, "_")
        if length(instName)==1
            print(fout, file * " & ")
        else
            ss = ""
            for ele in instName[1:end-1]
                if ele != "" ss = ss * ele * "\\_" end
            end
            if instName[end] != "" ss = ss * instName[end] end 
            print(fout, ss * " & ")
        end
        
        times = []
        pts = [] ; X = []

        # write dichotomy result 
        include(work_dir * "/dicho/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")
        push!(times, total_times_used); push!(pts, size_Y_N)

        # write ϵ-constraint result (ϵ = 0.5 by default)
        if isfile(work_dir * "/epsilon/" * file)
            include(work_dir * "/epsilon/" * file)
            push!(times, total_times_used); push!(pts, size_Y_N)
        else
            push!(times, -1); push!(pts, -1)
        end

        # write B&B result 
        if isfile(work_dir * "/bb/" * file)
            include(work_dir * "/bb/" * file)
            push!(times, total_times_used); push!(pts, size_Y_N) 
            push!(X, size_X_E)
        else
            push!(times, -1); push!(pts, -1)
            push!(X, -1)
        end

        # write B&C result 
        if isfile(work_dir * "/bc/" * file)
            include(work_dir * "/bc/" * file)
            push!(times, total_times_used); push!(pts, size_Y_N) 
            push!(X, size_X_E)
        else
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
    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "\\label{tab:table_compare_$instances }")
    println(fout, "\\end{table}")
    close(fout)

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
    avg_dicho = 0.0 ; avg_dominance = 0.0 ; avg_incumbent = 0.0
    avg_totalNodes = 0.0 ; avg_pruned = 0.0 ; n=0

    for file in readdir(dir * "/bb/")
        if split(file, ".")[end] == "png"
            continue
        end
        n += 1

        instName = split(file, "_")
        if length(instName)==1
            print(fout, file * " & ")
        else
            ss = ""
            for ele in instName[1:end-1]
                if ele != "" ss = ss * ele * "\\_" end
            end
            if instName[end] != "" ss = ss * instName[end] end 
            print(fout, ss * " & ")
        end

        include(dir * "/bb/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")
        print(fout, string(total_times_used)* " & " * string(relaxation_time) * " & " *
            string(test_dominance_time) * " & " * string(update_incumbent_time) * " & " *
            string(total_nodes) * " & " * string( pruned_nodes) * " & " * string(tree_size) * " & " *
            string(size_Y_N) * " & " * string(size_X_E)
        )

        avg_dicho += relaxation_time ; avg_dominance += test_dominance_time ; avg_incumbent += update_incumbent_time
        avg_totalNodes += total_nodes ; avg_pruned += pruned_nodes

        println(fout, "\\\\")
    end

    avg_dicho = round(avg_dicho/n, digits = 2) ; avg_dominance = round(avg_dominance/n, digits = 2) ; avg_incumbent = round(avg_incumbent/n, digits = 2) ;
    avg_totalNodes = round((avg_totalNodes-avg_pruned)/n, digits = 2) ; avg_pruned = round(avg_pruned/n, digits = 2) ;

    latex = raw"""\bottomrule
    \end{tabular}%
    }%
    \caption{The detailed experimental information about BO01B\&B algorithm.}
    \label{tab:table_bb}
    \end{table}
    """
    println(fout, latex)
    close(fout)

    # # Creating plot
    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie([avg_dicho, avg_dominance, avg_incumbent],
            explode=[0.1, 0.2, 0.3], labels = ["BOLP", "dominance", "incumbent"],
            autopct="%1.1f%%", shadow=false, startangle=45
        )

    # Adding legend
    ax.legend(wedges, ["BOLP", "dominance", "incumbent"],
            title ="Components",
            loc ="upper left",
            fontsize=7)
    
    plt.setp(autotexts, size = 8, weight ="bold")
    ax.set_title("CPU time of each components")
    
    savefig(dir * "/pieChartB&BTimes.png")
    plt.close()

    # # Creating plot
    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie([avg_totalNodes, avg_pruned],
            explode=[0.1, 0.2], labels = ["not fathomed", "pruned"],
            autopct="%1.1f%%", shadow=true, startangle=60
        )

    # Adding legend
    ax.legend(wedges, ["not fathomed", "pruned"],
            title ="#Nodes",
            loc ="upper right",
            fontsize=7)
    
    plt.setp(autotexts, size = 8, weight ="bold")
    ax.set_title("Number of nodes")
    
    savefig(dir * "/pieChartB&BNodes.png")
    plt.close()
end


function MOBC_perform(instances::String)
    bc = "/bc"
    dir = "../../results/" * instances * "/" * bc
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
    avg_dicho = 0.0 ; avg_pool = 0.0 ; avg_sepa = 0.0 ; avg_cuts = 0.0
    avg_totalNodes = 0.0 ; avg_pruned = 0.0 ; n=0
    avg_sp = 0.0 ; avg_mp = 0.0

    for file in readdir(dir * "/")
        if split(file, ".")[end] == "png"
            continue
        end
        n += 1

        instName = split(file, "_")
        if length(instName)==1
            print(fout, file * " & ")
        else
            ss = ""
            for ele in instName[1:end-1]
                if ele != "" ss = ss * ele * "\\_" end
            end
            if instName[end] != "" ss = ss * instName[end] end 
            print(fout, ss * " & ")
        end

        include(dir * "/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")
        print(fout, string(total_nodes) * " & " * string(pruned_nodes) * " & " * string(ite_total) * " & "*
            string( round(ite_total/total_nodes, digits = 2) ) * " & " * string(cuts_applied) * " & " * string(sp_cuts)*
            " & " * string(mp_cuts) * " & " * string(cuts_total) * " & " * string(times_total_for_cuts) * " & " *
            string(times_calling_dicho) * " & " * string(times_oper_cutPool) * " & " *string(times_calling_separators) * " & "*
            string(times_add_retrieve_cuts) * " & " *string(total_times_used) * " & " *string(size_Y_N)
        )

        avg_dicho += times_calling_dicho ; avg_pool += times_oper_cutPool ; avg_sepa += times_calling_separators ; avg_cuts += times_add_retrieve_cuts
        avg_totalNodes += total_nodes ; avg_pruned += pruned_nodes
        avg_sp += sp_cuts ; avg_mp += mp_cuts

        println(fout, "\\\\")
    end

    avg_dicho = round(avg_dicho/n, digits = 2) ; avg_pool = round(avg_pool/n, digits = 2) ; avg_sepa = round(avg_sepa/n, digits = 2) ;avg_cuts = round(avg_cuts/n, digits = 2) ;
    avg_totalNodes = round((avg_totalNodes-avg_pruned)/n, digits = 2) ; avg_pruned = round(avg_pruned/n, digits = 2) ;
    avg_sp = round(avg_sp/n, digits = 2) ; avg_mp = round(avg_mp/n, digits = 2) ;


    latex = raw"""\bottomrule
    \end{tabular}%
    }%
    \caption{.}
    \label{tab:table_bc}
    \end{sidewaystable}
    """
    println(fout, latex)
    close(fout)

    # # Creating plot
    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie([avg_dicho, avg_pool, avg_sepa, avg_cuts],
            explode=[0.1, 0.2, 0.3, 0.4], labels = ["BOLP", "cutpool", "separator", "cuts oper"],
            autopct="%1.1f%%", shadow=false, startangle=45
        )

    # Adding legend
    ax.legend(wedges, ["BOLP", "cutpool", "separator", "cuts oper"],
            title ="Components",
            loc ="upper left",
            fontsize=7)
    
    plt.setp(autotexts, size = 8, weight ="bold")
    ax.set_title("CPU time of each components")
    
    savefig(dir * "/pieChartB&CTimes.png")
    plt.close()

    # # Creating plot
    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie([avg_totalNodes, avg_pruned],
            explode=[0.1, 0.2], labels = ["not fathomed", "pruned"],
            autopct="%1.1f%%", shadow=true, startangle=60
        )

    # Adding legend
    ax.legend(wedges, ["not fathomed", "pruned"],
            title ="#Nodes",
            loc ="upper right",
            fontsize=7)
    
    plt.setp(autotexts, size = 8, weight ="bold")
    ax.set_title("Number of nodes")
    
    savefig(dir * "/pieChartB&CNodes.png")
    plt.close()

    # # Creating plot
    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie([avg_sp, avg_mp],
            explode=[0.2, 0.1], labels = ["single-point", "multi-point"],
            autopct="%1.1f%%", shadow=true, startangle=120
        )

    # Adding legend
    ax.legend(wedges, ["single-point", "multi-point"],
            title ="#Cuts",
            loc ="upper right",
            fontsize=7)
    
    plt.setp(autotexts, size = 8, weight ="bold")
    ax.set_title("Number of cuts")
    
    savefig(dir * "/pieChartB&CCuts.png")
    plt.close()

end


comparisonMethods("momhMKPstu/MOBKP/set2")
detailedMOBB_perform("momhMKPstu/MOBKP/set2")
MOBC_perform("momhMKPstu/MOBKP/set2")