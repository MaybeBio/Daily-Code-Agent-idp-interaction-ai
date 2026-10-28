function GenerateSequenceByLine(Sequence; N=60)
    str = ""
    cnt = 1
    M= ceil(Int32,length(Sequence)/N)
    for I in 1:M
        str *="$(1+(I-1)*N) & \\texttt{"
        for j in 1:N÷10
            for i in 1:10
                str *= Sequence[cnt]
                cnt += 1
                if cnt > length(Sequence) 
                    str *= "}\\\\ \n "
                    return str,  M+2
                end
            end
            str *= " "
        end
        str *= "}\\\\ \n "
    end
    return str, M+2
end

function GenerateSubtableWrapper(text; Caption="")
    return "\\vspace{0.1cm}
    \\begin{subtable}{\\textwidth}
    \\caption{$Caption}
    %\\centering
    \\begin{tabular}{r|l}
    $text    
    \\end{tabular}
    \\end{subtable}"
end

function GenerateSequenceTables(Sequences, Names, Outpath;  LineMax=46)
    TableContents = Vector{String}([])
    LinesPerTable=Vector{Int32}([])
    NSeq= length(Sequences)
    for (cnt,Seq) in enumerate(Sequences)
        cont, lines= GenerateSequenceByLine(Seq)
        push!(TableContents, cont)
        append!(LinesPerTable, lines)
    end

    Pages = Vector{Vector{String}}([[]])
    cnt =0
    page_cnt=1
    for i in 1:NSeq
        if (cnt+ LinesPerTable[i])<=LineMax
            push!(Pages[page_cnt], TableContents[i])
            cnt +=  LinesPerTable[i]
        else
            cnt = 0 #LinesPerTable[i]
            page_cnt += 1
            push!(Pages,[TableContents[i]])
        end
    end

    cnt = 1
    for (PID, page) in enumerate(Pages)
        out ="\\begin{table}[htb]\n
        \\centering"
        for table in page
            out *= GenerateSubtableWrapper(table; Caption="\\textbf{$(Names[cnt]):}")
            cnt += 1
        end
        out *=  "\\label{tab:App_SEQ_$(PID)}\n \\end{table}"
        write("$Outpath/Sequence_Table_$PID.tex", out)
    end
    return nothing
end