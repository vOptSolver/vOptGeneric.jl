

function comparisonThreeMethods(instances::String)
    dir = "../../results/" * instances
    @assert isdir(dir) "This directory doesn't exist $dir !"

    fout = open(dir * "/comparisonTable.tex", "w")

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

    for file in readdir(dir * "/dicho")
        if split(file, ".")[end] == "png"
            continue
        end

        print(fout, file * " & ")
        times = []
        pts = []

        # write dichotomy result 
        include(dir * "/dicho/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")
        # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
        push!(times, total_times_used); push!(pts, size_Y_N)

        # write ϵ-constraint result (ϵ = 1.0 by default)
        if isfile(dir * "/epsilon/epsilon_0.1/" * file)
            include(dir * "/epsilon/epsilon_0.1/" * file)
            # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
            push!(times, total_times_used); push!(pts, size_Y_N)
        else
            # print(fout, "- & - & ")
            push!(times, -1); push!(pts, -1)
        end

        # write B&B result 
        if isfile(dir * "/bb/default/" * file)
            include(dir * "/bb/default/" * file)
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
    % }%
    \caption{Comparison of the different algorithms performances.}
    \label{tab:table_compare}
    \end{table}
    """
    println(fout, latex)
    close(fout)
end




function detailedComparisonMethods(instances::String)
    dir = "../../results/" * instances
    @assert isdir(dir) "This directory doesn't exist $dir !"

    fout = open(dir * "/detailedComparisonMethods.tex", "w")

    latex = raw"""\begin{sidewaystable}[!h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m}  & \multicolumn{2}{c}{\textbf{Dichotomy}} & \multicolumn{9}{c}{\textbf{$\mathbf{\epsilon}$-constraint}} & \multicolumn{6}{c}{\textbf{Branch-and-bound}}
    \\
    \cmidrule(r){4-5} \cmidrule(r){6-14} \cmidrule(r){15-20}
    ~ & ~ & ~ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & $\mathbf{\epsilon}$ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & $\mathbf{\epsilon}$ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & $\mathbf{\epsilon}$ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} &  
    \textbf{fathom} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} &
    \textbf{fathom} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$}
    \\
    \midrule
    """
    println(fout, latex)

    for file in readdir(dir * "/dicho")
        if split(file, ".")[end] == "png"
            continue
        end

        print(fout, file * " & ")
        times = []
        pts = []

        # write dichotomy result 
        include(dir * "/dicho/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")
        print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")

        # write ϵ-constraint result 
        for ϵ in [0.1, 1.0, 5.0]
            print(fout, string(ϵ) * " & ")
            if isfile(dir * "/epsilon/epsilon_$ϵ/" * file)
                include(dir * "/epsilon/epsilon_$ϵ/" * file)
                print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
            else
                print(fout, "- & - & ")
            end
        end

        # write B&B result 
        print(fout, "yes & ")
        if isfile(dir * "/bb/default/" * file)
            include(dir * "/bb/default/" * file)
            print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
        else
            print(fout, "- & - & ")
        end
        print(fout, "no & ")
        if isfile(dir * "/bb/fathom_false/" * file)
            include(dir * "/bb/fathom_false/" * file)
            print(fout, string(total_times_used)* " & " * string(size_Y_N))
        else
            print(fout, "- & - ")
        end

        println(fout, "\\\\")
    end

    latex = raw"""\bottomrule
    \end{tabular}%
    }%
    \caption{The performance comparison between different algorithms with different settings.}
    \label{tab:table_compare2}
    \end{sidewaystable}
    """
    println(fout, latex)
    close(fout)
end


function detailedMOBB(instances::String)
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

    for file in readdir(dir * "/bb/default/")
        if split(file, ".")[end] == "png"
            continue
        end

        print(fout, file * " & ")
        times = []
        pts = []


        include(dir * "/bb/default/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")
        print(fout, string(total_times_used)* " & " * string(relaxation_time) * " & " *
            string(test_dominance_time) * " & " * string(update_incumbent_time) * " & " *
            string(total_nodes) * " & " * string( pruned_nodes) * " & " * string(tree_size) * " & " * string(size_Y_N)
        )

        println(fout, "\\\\")
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

detailedMOBB("smallExamples")
comparisonThreeMethods("smallExamples")
detailedComparisonMethods("smallExamples")