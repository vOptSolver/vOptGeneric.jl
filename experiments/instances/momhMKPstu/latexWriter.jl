
function comparisonThreeMethods(instances::String)
    dir = "../../results/" * instances
    @assert isdir(dir) "This directory doesn't exist $dir !"


    for class in readdir(dir)
        @assert isdir(dir * "/" * class) "This directory doesn't exist $dir/$class !"
        work_dir = dir * "/" * class

        fout = open(work_dir * "/comparisonTable.tex", "w")

        latex = raw"""\begin{table}[!h]
        \centering
        % \resizebox{\columnwidth}{!}{%
        \hspace*{-1cm}\begin{tabular}{lcccccccc}
        \toprule
        \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{Dichotomy}} & \multicolumn{2}{c}{\textbf{$\mathbf{\epsilon}$-constraint}}  & \multicolumn{2}{c}{\textbf{Branch-and-bound}}
        \\
        \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-9}
        ~ & ~ & ~ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} \\
        \midrule
        """
        println(fout, latex)
    
        for file in readdir(work_dir * "/dicho")
            if split(file, ".")[end] == "png"
                continue
            end
    
            print(fout, file * " & ")
            times = []
            pts = []
    
            # write dichotomy result 
            include(work_dir * "/dicho/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")
            # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
            push!(times, total_times_used); push!(pts, size_Y_N)
    
            # write ϵ-constraint result (ϵ = 0.5 by default)
            if isfile(dir * "/epsilon/" * file)
                include(dir * "/epsilon/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N)
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
            end
    
            # write B&B result 
            if isfile(dir * "/bb/" * file)
                include(dir * "/bb/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N)
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
            end
    
            for i=1:3
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] == minimum(times)
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    i<3 ? print(fout, " - & ") : print(fout, " - ")
                elseif pts[i] == maximum(pts)
                    i < 3 ? print(fout, " \\textbf{" * string(pts[i]) * "} & ") : print(fout, " \\textbf{" * string(pts[i]) * "} ")
                else
                    i<3 ? print(fout, string(pts[i]) * " & ") : print(fout, string(pts[i]) * " ")
                end
            end
    
            println(fout, "\\\\")
        end
    
        latex = raw"""\bottomrule
        \end{tabular}
        % }%"""
        println(fout, latex)
        println(fout, "\\caption{Comparison of the different algorithms performances for instances $class .}")
        println(fout, "\\label{tab:table_compare_$class }")
        println(fout, "\\end{table}")
        close(fout)
    end

end


function detailedMOBB_perform(instances::String)
    dir = "../../results/" * instances
    @assert isdir(dir) "This directory doesn't exist $dir !"

    fout = open(dir * "/detailedMOBB.tex", "w")

    latex = raw"""\begin{table}[!h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{4}{c}{\textbf{Time(s)}} & \multicolumn{2}{c}{\textbf{Nodes}}  & \textbf{Tree(MB)} & \textbf{$|\mathcal{Y}_N|$}
    \\
    \cmidrule(r){4-7} \cmidrule(r){8-9} 
    ~ & ~ & ~ & \textbf{total} &\textbf{relax} & \textbf{dominance} & \textbf{incumbent} & \textbf{total} & \textbf{pruned} & ~ & ~\\
    \midrule
    """
    println(fout, latex)

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
                string(total_nodes) * " & " * string( pruned_nodes) * " & " * string(tree_size) * " & " * string(size_Y_N)
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

        println(fout, "avg & " * string(avg_n) * " & " * string(avg_m) * " & " * string(avg_totalT) * "  &  " *
            string(avg_relaxT) * " & " * string(avg_dominanceT) * " & " * string(avg_incumbentT) * " & " *
            string(avg_totalN) * " & " * string(avg_prunedN) * " & " * string(avg_treeSize) * " & " * string(avg_YN) * " \\\\ \\cline{1-11}")
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
end

detailedMOBB_perform("momhMKPstu/MOBKP")



# comparisonThreeMethods("momhMKPstu")