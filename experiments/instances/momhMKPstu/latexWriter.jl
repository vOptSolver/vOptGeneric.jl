
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

comparisonThreeMethods("momhMKPstu")