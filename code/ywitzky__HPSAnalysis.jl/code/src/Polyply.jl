
module Polyply
using HPSAnalysis
using ..BioData: OneToThree
using Scratch

polyply = "$(HPSAnalysis.EnvironmentPath)/bin/polyply"

function ConvertGroToPDB(Path, Filename)
    GMX_Path="/localscratch/Programs/gromacs-2022.2/bin/gmx"
    run(`$GMX_Path editconf -f $(Path)$(Filename).gro -o $(Path)$(Filename).pdb`)
end

function GenerateITPFilesOfSequence(Names, Sequences, OutputPath)
    NameSet = Set(Names)
    for name in (NameSet)
        ind = findfirst(x-> x==name, Names)
        SeqString = []

        for AA in Sequences[ind]
            ### treat phosphorylated resdidues like 
            if AA == '#'
                push!(SeqString, "$(OneToThree['S']):1")
            elseif AA=='&'
                push!(SeqString, "$(OneToThree['T']):1")
            elseif AA=='@'
                push!(SeqString, "$(OneToThree['Y']):1")
            else
                push!(SeqString, "$(OneToThree[AA]):1")
            end
        end
        run(`$polyply gen_params  -name $(name) -seq $SeqString -lib martini3  -o $(OutputPath)$(name).itp`) 
    end
end

function GenerateSlabTopologyFile(Filename, ITPPath, Names, SimulationName)
    f = open(Filename,"w+")

    write(f,"#include \"$(HPSAnalysis.PkgPath)/data/Polyply/martini_v3.0.0.itp\" \n")

    write(f, "[ atomtypes ]\n")
    write(f, "VS 0.00 0.000 V 0.0 0.0\n\n")


    write(f, "[ nonbond_params ]\n")
    write(f, "VS    W   1 0.4650000000    0.5000000000\n\n")
    
    NameSet = Set(Names)
    Occurences = counter(Names)
    for name in NameSet
        write(f, "#include \"$(ITPPath)$(name).itp\" \n")
    end
    write(f,"\n[ system ]\n; name\nProtein Sim\n\n[ molecules ]\n; name  number\n")

    for (ind, name) in enumerate(NameSet)
        write(f, "$name $(Occurences[name])\n")
    end
    close(f)
end

function GenerateCoordinates(SimulationPath, SimulationName, Box, TopologyFile)
    run(`$polyply gen_coords -p $TopologyFile -o $(SimulationPath)$(SimulationName).gro -name $(SimulationName) -box $(Box)`)
end

function readSimpleGRO(Filename, x,y,z)
    f = open(Filename)
    readline(f) ### command that generate the file
    readline(f) ### the number of lines
    ind=1
    while !eof(f) 
        line = readline(f)
        if length(line)>15 && line[14:15]=="CA"
            x[ind, 1] = parse(Float32, line[21:28])
            y[ind, 1] = parse(Float32, line[29:37])
            z[ind, 1] = parse(Float32, line[38:44])
            ind+=1
        end
    end
    ### convert from nm to Angstroem
    x .*= 10.0
    y .*= 10.0
    z .*= 10.0
end

function readPDB(Filename, x,y,z)
    f = open(Filename)
    ind =1
    while !eof(f) 
        line = readline(f)
        if length(line)<6
            continue
        elseif line[1:6]=="ATOM  " 
            if line[14:15]=="BB"
                x[ind, 1] = parse(Float32, line[31:38])
                y[ind, 1] = parse(Float32, line[39:46])
                z[ind, 1] = parse(Float32, line[47:54])
                ind += 1
            end
        end
    end
end

end 