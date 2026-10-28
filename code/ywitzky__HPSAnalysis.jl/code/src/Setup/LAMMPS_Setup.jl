@doc raw"""
    writeHPSLammps(fileName,Sequences,AtomTypes,LongAtomTypes,LongAtomTypesToRes,InitStyle,ChargeTemperSteps,ChargeTemperSwapSteps,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,pH,NAtoms,NBonds,NAngles,NDihedrals,Info,NAtomTypes)

Writes the datas for the simulation in different files.
    
**Arguments**
- from writeStartConfiguration()

**Creates**:
* Writes the datas for the simulation.
"""
function writeHPSLammps(fileName,Sequences,AtomTypes,LongAtomTypes,LongAtomTypesToRes,InitStyle,ChargeTemperSteps,ChargeTemperSwapSteps,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,pH,NAtoms,NBonds,NAngles,NDihedrals,Info,NAtomTypes)
    writeHPSLammpsScript( fileName*".lmp",StartFileName, AtomTypes, LongAtomTypes, AaToId, LongAtomTypesToRes, OneToCharge, OneToSigma, OneToLambda, dihedral_eps, InitStyle, SimulationType, Temperature, AlphaAddition, false, NSteps; SaltConcentration=SaltConcentration, pH=pH, ChargeTemperSteps=ChargeTemperSteps, ChargeTemperSwapSteps=ChargeTemperSwapSteps,WriteOutFreq=WriteOutFreq)
    writeHPSLammpsScript( fileName*"_restart.lmp",StartFileName, AtomTypes, LongAtomTypes, AaToId, LongAtomTypesToRes, OneToCharge, OneToSigma, OneToLambda, dihedral_eps, InitStyle, SimulationType, Temperature, AlphaAddition, true, NSteps; SaltConcentration=SaltConcentration, pH=pH, ChargeTemperSteps=ChargeTemperSteps, ChargeTemperSwapSteps=ChargeTemperSwapSteps,WriteOutFreq=WriteOutFreq)

    file = open(StartFileName*".txt", "w");
    write(file, Info)

    write(file, "\n\t $NAtoms \t atoms\n")
    write(file, "\t $NBonds \t bonds\n")
    write(file, "\t $NAngles \t angles\n")
    write(file, "\t $NDihedrals \t dihedrals\n\n")

    write(file, "\t $NAtomTypes \t atom types\n")
    write(file, "\t 1 \t bond types\n")
    write(file, "\t 1 \t angle types\n")
    write(file, "\t $NDihedralsTypes \t dihedral types\n\n")

    write(file,"\t $(BoxSize[1]) \t $(BoxSize[2]) \t xlo xhi\n")
    write(file,"\t $(BoxSize[3]) \t $(BoxSize[4]) \t ylo yhi\n")
    write(file,"\t $(BoxSize[5]) \t $(BoxSize[6]) \t zlo zhi\n\n")

    write(file, "Masses \n#\n")
    for (index, value) in enumerate(LongAtomTypes)
        write(file, "\t $(AaToId[value]) \t $(OneToMass[value]) \t ### $(value in AtomTypes ? value : LongAtomTypesToRes[value] ) \n")
    end

    write(file,"\nAtoms\n# A comment that is needed to read stuff accurately\n")

    atomid=0
    moleculeid=0
    for (SeqId, Sequence) in enumerate(Sequences)
        moleculeid+=1
        for (ResId,Res) in enumerate(Sequence)
            atomid +=1
            write(file, "\t $atomid \t $(@sprintf("%6i", moleculeid)) \t  $(@sprintf("%2i", AaToId[Res])) \t  $(@sprintf("%3.5f", OneToCharge[Res])) \t $(@sprintf("%5.5f", pos[SeqId,ResId,1])) \t $(@sprintf("%5.5f", pos[SeqId,ResId,2])) \t $(@sprintf("%5.5f", pos[SeqId,ResId,3])) \t $(@sprintf("%i", image[SeqId,ResId,1])) \t $(@sprintf("%i", image[SeqId,ResId,2])) \t $(@sprintf("%i", image[SeqId,ResId,3]))\n")
        end
    end 

    write(file,"\nBonds\n#\n")
    bonds = getBonds(Sequences;M=2)
    for bid in axes(bonds,1)
        write(file, "\t $bid \t 1 \t  $(bonds[bid, 1]) \t  $(bonds[bid, 1]) \n")
    end
    write(file, "\n\n")


    write(file,"\nAngles\n#\n")

    angles = getBonds(Sequences;M=3)
    for bid in axes(bonds,1)
        write(file, "\t $bid \t 1 \t  $(angles[bid, 1]) \t  $(angles[bid, 1]) \t  $(angles[bid, 2])\n")
    end
    write(file, "\n\n")


    write(file,"\nDihedrals\n#\n")
    atomid=0
    dihedralid = 0
    for (SeqId, Seq) in enumerate(Sequences)
        for (ResId,Res) in enumerate(Seq)
            atomid +=1
            if (ResId>(length(Seq)-3) )
                continue
            else
                dihedralid += 1
                if MixingRule=="1-1001-1"
                    Res_min = (ResId-1)>=1 ? AaToId[Seq[ResId-1]] : 0
                    Res_max = (ResId)<=(length(Seq)-4) ? AaToId[Seq[ResId+4]] : -1
                    key = (sort([Res_min,AaToId[Res], AaToId[Seq[ResId+3]], Res_max]))
                elseif MixingRule=="1001"
                    key = (sort([AaToId[Res], AaToId[Seq[ResId+3]]]))
                elseif MixingRule=="0110"
                    key = (sort([AaToId[Seq[ResId+1]], AaToId[Seq[ResId+2]]]))
                end
                write(file, "\t $dihedralid \t $(dihedral_long_map[key]) \t  $atomid \t  $(atomid+1) \t $(atomid+2) \t $(atomid+3)\n")
            end
        end
    end
    close(file)    
end


function writeSlurmScript(fileName, Proteins, Path, RunsPerProtein, OMPCores, Restart=false)
    file = open(fileName, "w");
    write(file,"#!/bin/bash\nexport lmp=/p/home/jusers/witzky1/juwels/LAMMPS_HREMD/bin/lmp\n\n")
    for (ProtId, Protein) in enumerate(Proteins)
        for Run in 1:RunsPerProtein[ProtId]
            write(file,"cd ./$(Protein)/RUN_$( lpad(Run,3,"0"))/\n")
            if (Restart)
                write(file,"srun -n 1 --exclusive --exact --cpus-per-task=$(OMPCores[ProtId]) \$lmp -in ./$(Protein)_restart.lmp -sf omp -package omp $(OMPCores[ProtId]) &\n")
            else
                write(file,"srun -n 1 --exclusive --exact --cpus-per-task=$(OMPCores[ProtId]) \$lmp -in ./$(Protein).lmp -sf omp -package omp $(OMPCores[ProtId]) & \n")
            end
            write(file, "cd ../../\n\n")
        end
    end
    write(file, "wait ")
    close(file)
end

function writeCollectedSlurmScript(Path, Proteins, RelPaths,MPICores,OMPCores; ProteinsPerSlurmFile=2, Restarts=6, SlurmAccName="rsproteins", LmpAddOn="", CoresPerNodes=48, JobName="RS")
    Nodes=Int32(ceil(maximum(MPICores)*ProteinsPerSlurmFile/CoresPerNodes))
    slurm_file = open(Path*"run_all.sh", "w");
    write(slurm_file,"#!/bin/bash\nACCOUNT=$(SlurmAccName)\nRUNTIME=\"24:00:00\"\nPARTITION=\"batch\"\nNAME=\"$JobName\"\nMAIL=\" --mail-type=END,FAIL --mail-user=ywitzky@students.uni-mainz.de\"\nmodule load Intel\nmodule load IntelMPI\nmodule load HDF5/1.14.2-serial\nexport BASEPATH=\$PWD\n\n\n")
    cnt=0
    StartFile = open(Path*"test.sh", "w")
    RestartFile = open(Path*"test.sh", "w")
    for (RunNum, Protein) in enumerate(Proteins)
        if(mod(RunNum,ProteinsPerSlurmFile)==1  || ProteinsPerSlurmFile==1)
            if(RunNum÷ProteinsPerSlurmFile>0)
                write(StartFile, "\n\nwait ")
                write(RestartFile, "\n\nwait ")
                close(StartFile)
                close(RestartFile)
                if ProteinsPerSlurmFile>1
                    write(slurm_file,"\njid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL StartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
                    for ind in 1:Restarts
                        write(slurm_file,"jid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL --dependency=afterok:\${jid$(cnt-1)##* } RestartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
                    end
                end
            end
            StartFile = open(Path*"StartScript_$(RunNum÷ProteinsPerSlurmFile)", "w");
            write(StartFile,"#!/bin/bash\nexport lmp=/p/home/jusers/witzky1/juwels/LAMMPS_HREMD/bin/lmp\n\n")
            RestartFile = open(Path*"RestartScript_$(RunNum÷ProteinsPerSlurmFile)", "w");
            write(RestartFile,"#!/bin/bash\nexport lmp=/p/home/jusers/witzky1/juwels/LAMMPS_HREMD/bin/lmp\n\n")
        end
        write(StartFile,"cd ../$(RelPaths[RunNum]) \n")
        write(RestartFile,"cd ../$(RelPaths[RunNum]) \n")\
        write(RestartFile,"srun -n $(MPICores[RunNum]) --exclusive --exact --cpus-per-task=$(OMPCores[RunNum]) \$lmp -in \$BASEPATH/../$(RelPaths[RunNum])/$(Protein)_restart.lmp $(LmpAddOn)  -sf omp -package omp $(OMPCores[RunNum]) &\n")
        write(StartFile,"srun -n $(MPICores[RunNum]) --exclusive --exact --cpus-per-task=$(OMPCores[RunNum]) \$lmp -in \$BASEPATH/../$(RelPaths[RunNum])/$(Protein).lmp $(LmpAddOn)  -sf omp -package omp $(OMPCores[RunNum]) & \n")
        write(StartFile, "cd \$BASEPATH\n\n")
        write(RestartFile, "cd \$BASEPATH\n\n")
    end

    if(length(Proteins)÷ProteinsPerSlurmFile>0 || ProteinsPerSlurmFile==1 )

        RunNum=length(Proteins)+1
        write(StartFile, "\n\nwait ")
        write(RestartFile, "\n\nwait ")
        close(StartFile)
        close(RestartFile)
        write(slurm_file,"\njid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL StartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
        for ind in 1:Restarts
            write(slurm_file,"jid$(cnt+=1)=\$(sbatch -J Calva2_$(RunNum÷ProteinsPerSlurmFile-1) -A \$ACCOUNT -p \$PARTITION -t \$RUNTIME \$MAIL --dependency=afterok:\${jid$(cnt-1)##* } RestartScript_$(RunNum÷ProteinsPerSlurmFile-1))\n")
        end
        #end
    end
    
    close(slurm_file)
end

function writeHPSLammpsScript(fileName, StartFileName, AtomTypes, LongAtomTypes, AaToId, LongAtomTypesToRes, OneToCharge, OneToSigma, OneToLambda, Dihedral_eps, InitStyle, SimulationType, Temperature, AlphaStruct, Restart, NTimeSteps; WriteOutFreq=100_000,SaltConcentration=-1, pH = 6, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, OutFormat="h5md")

    ChargeTemperSim=length(ChargeTemperSteps)>0

    (ϵ_r, κ) = DetermineYukawaInteractions(;SimulationType=SimulationType, Temperature=Temperature, SaltConcentration=SaltConcentration)

    file = open(fileName, "w")
    write(file,"
    ### Initialisation

    units real ### https://docs.lammps.org/stable/units.html
    dimension 3
    newton on ### might be faster if off depends on parallelization https://docs.lammps.org/stable/newton.html
    processors * * * grid numa #map zyx
    ### processors pre dimension, can optimize a lot here.
    boundary p p p # periodic in each dimension 
    atom_style full 

    atom_modify id yes  #ids will be assigned to each atom 
    #atom_modify map hash 
    #atom_modify 10000 2 #resort ids every 10000 steps for cache improvements 

    # Reading file i easier\n")


    if ChargeTemperSim
        write(file,"\n    ### set up H-REMD
    variable q world $(join(string.(ChargeTemperSteps).*" "))
    variable id world $(join(string.(collect(1:length(ChargeTemperSteps))).*" "))
    variable energyfile string \"energy_\${id}.dat\" \n")
    else 
        write(file," variable energyfile string \"energy.dat\" \n")
    end

    if (Restart)
    write(file, "read_restart ./sim$(ChargeTemperSim ?  "_\${id}" : "").restart \n")
    else
    write(file, "read_data $(StartFileName)$(ChargeTemperSim ?  "_\${id}" : "").txt\n")
    end

    write(file, "neighbor 3.5 multi\n")
    if SimulationType=="HPS-Alpha"
        write(file,"neigh_modify every 5 delay 0
    neigh_modify collection/interval 5 20 22 24 26 36\n")
    elseif SimulationType=="Calvados2"
        write(file,"neigh_modify every 5 delay 0
    neigh_modify collection/interval 2 20.5 40.5\n")
    end

    write(file,"
    bond_style  harmonic")
    if AlphaStruct
        write(file,"
    angle_style bch
    dihedral_style gaussian")
    end
    write(file,"

    dielectric  $(@sprintf("%.5f", ϵ_r))
    
    # bond potential parameters
    bond_coeff          1   $(@sprintf("%.5f", get!(BondStrength, SimulationType,10.000)))    3.800000") 

    if AlphaStruct
        write(file,"

        # angle potential parameters 
        angle_coeff         1    4.300000

        # dihedral angle parameters\n")
        for (i, eps) in enumerate(Dihedral_eps)
            write(file,"dihedral_coeff        $(i)    $(eps)\n")
        end
    end
    write(file,"\n
    ### Coulomb and Ashbaugh Hatch Potential
    pair_style  ljlambda $(@sprintf("%.5f", 1/κ)) 0.0 $(Cutoff[SimulationType])
    dielectric  $(@sprintf("%.5f", ϵ_r))\n")
    if AlphaStruct
        write(file,"special_bonds lj/coul 0.0 0.0 0.0\n")
    else
        write(file,"special_bonds lj/coul 0.0 1.0 1.0\n")
    end
    #    special_bonds ljlambda 0.0 1.0 1.0\n")

    #for (id, ResType) in enumerate(AtomTypes)
    #    write(file, "pair_coeff $id $id $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType])  $(BioData.OneToHPSCharge[ResType])\n")
    #end



    for  ResType in LongAtomTypes
        id = AaToId[ResType]
        for  ResType2 in LongAtomTypes
            id2 = AaToId[ResType2]

            if id2<id
                continue
            end
            charge = OneToCharge[ResType]*OneToCharge[ResType2]
            lambda = (OneToLambda[ResType]+OneToLambda[ResType2])/2.
            sigma = (OneToSigma[ResType]+OneToSigma[ResType2])/2.
            if SimulationType=="HPS-Alpha" ||  SimulationType=="ArashModell"
                cutoff=4*sigma
                if abs(charge) >0.001 
                    cutoff2= Cutoff[SimulationType]
                else
                    cutoff2= 0
                end
            elseif SimulationType=="Calvados2"
                cutoff=20
                if abs(charge) >10^-8
                    cutoff2= 40
                else
                    cutoff2= 0
                end
            end
            write(file, "pair_coeff $id $id2  $(@sprintf("%.5f", get!(EpsilonAshbaughHatch, SimulationType,0.2000))) $(@sprintf("%.5f", sigma)) $(@sprintf("%.5f", lambda)) $(cutoff) $(cutoff2)\n")
        end
        #write(file, "pair_coeff $id $id ah $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType]) \n")
    end
    #write(file, "\npair_modify mix arithmetic shift yes\n")
    Charged = ""
    Uncharged = ""
    for (id, ResType) in enumerate(LongAtomTypes)
        if OneToCharge[ResType]!=0
            Charged *= "$id "
        else
            Uncharged *= "$id "
        end
    end
    if length(Charged)>0 && length(Uncharged)>0
    write(file, "group charged type $Charged
        group uncharged type $Uncharged
        
        comm_style tiled
        fix 4 all balance 100000 1.1 rcb weight group 2 charged 8 uncharged 1
        balance 1.1 rcb weight group 2 charged 8 uncharged 1\n        ")
    end

    if InitStyle=="Slab"
        
        write(file, "#comm_style tiled
    #fix 4 all balance 50000 1.05 rcb\n")
    end

    
    if ~(Restart)
        write(file,"\nminimize 1.0e-4 1.0e-6 500 1000\n")
        write(file, "reset_timestep 0\n")
    end

    write(file, "run_style verlet\n")

    if OutFormat=="xyz"
        write(file, "dump trajectory all xyz $WriteOutFreq ./Trajectory$(ChargeTemperSim ?  "_\${id}" : "").xyz")
    elseif OutFormat=="xtc"
        write(file,"dump trajectory all xtc $WriteOutFreq ./Trajectory$(ChargeTemperSim ?  "_\${id}" : "").xtc\n")
    elseif OutFormat=="h5md"
        write(file,"dump trajectory all h5md $WriteOutFreq ./Trajectory$(ChargeTemperSim ?  "_\${id}" : "").h5 position image create_group yes\n")
    end
    if (Restart)
        write(file, "dump_modify trajectory append yes\n")
    end

    
    write(file,"
    fix 1 all nve 
    fix 2 all langevin $(Temperature) $(Temperature) $(SimulationType=="Calvados2" || SimulationType=="ArashModell" ? 100 : 1000) $(rand(1:100_000)) ### Tstart, Tstop, Dampening in fs  (1ps dampening in paper), seed\n")
    if (Restart)
    write(file,"
    fix 3 all print $WriteOutFreq  \"\$(step) \$(time) \$(temp) \$(press) \$(etotal) \$(pe) \$(ke) \$(evdwl) \$(ecoul) \$(epair) \$(ebond) \$(eangle) \$(edihed)\" append \${energyfile} screen no\n")
    else
        write(file,"
    fix 3 all print $WriteOutFreq  \"\$(step) \$(time) \$(temp) \$(press) \$(etotal) \$(pe) \$(ke) \$(evdwl) \$(ecoul) \$(epair) \$(ebond) \$(eangle) \$(edihed)\" file \${energyfile} screen no\n")
    end

    write(file,"   
    restart 100000000 ./Restart/BackUp$(ChargeTemperSim ?  "_\${id}" : "").*
    timestep 10.0
    timer timeout 23:30:00")

    if ChargeTemperSim
        write(file,"   
    chargetemper $NTimeSteps $ChargeTemperSwapSteps \$q 2 $(rand(1:100_000)) $(rand(1:100_000)) ### timesteps swapfreq chargeinit temp_fix_id seed1 seed2 ")
    else
    write(file,"   
    run $(NTimeSteps) upto")
    end
    write(file,"
    write_restart ./sim$(ChargeTemperSim ?  "_\${id}" : "").restart")
    close(file)
end

function writeLammpsScript(fileName, AtomTypes; WriteOutFreq=100000)
    file = open(fileName, "w")

    write(file," ### Initialisation

    units real ### https://docs.lammps.org/stable/units.html
    dimension 3
    newton on ### might be faster if off depends on parallelization https://docs.lammps.org/stable/newton.html
    processors * * * grid numa #map zyx
    ### processors pre dimension, can optimize a lot here.
    boundary p p p # periodic in each dimension 
    atom_style full 
    bond_style harmonic/omp

    atom_modify id yes  #ids will be assigned to each atom 
    #atom_modify map hash 
    #atom_modify 10000 2 #resort ids every 10000 steps for cache improvements 


    # Reading file i easier
    read_data ./Start_conf.txt
    neighbor 2.0 multi
    neigh_modify every 5
    neigh_modify collection/interval 5 20 22 24 26 36


    ### Simulation Settings
    #Bonds
    bond_coeff 1 10.0 3.82   ### ????, k=10 kcal/mol/AA^2, sigma_0=3.82AA


    ### Coulomb and Ashbaugh Hatch Potential
    pair_style ah/cut/coul/cut 35 0.1 0.2 ### Dielectric constant
    dielectric 80.0 		### dielectric constant \n")

    for (id, ResType) in enumerate(AtomTypes)
        write(file, "pair_coeff $id $id $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType])  $(BioData.OneToHPSCharge[ResType])\n")
    end

    #=for (id, ResType) in enumerate(AtomTypes)
        for (id2, ResType2) in enumerate(AtomTypes)
            if id2<id
                continue
            end
            write(file, "pair_coeff $id $id2 ah $((BioData.OneToHPSUrrySigma[ResType]+BioData.OneToHPSUrrySigma[ResType2])/2.) $((BioData.OneToHPSUrryLambda[ResType]+BioData.OneToHPSUrryLambda[ResType])/2.) \n")
        end
        #write(file, "pair_coeff $id $id ah $(BioData.OneToHPSUrrySigma[ResType]) $(BioData.OneToHPSUrryLambda[ResType]) \n")
    end=#
    write(file, "\npair_modify mix arithmetic shift yes\n")


    write(file, "run_style verlet

    angle_style none
    dihedral_style none
    dump trajectory2 all xyz $(WriteOutFreq) ./Trajectory.xyz 


    minimize 1.0e-4 1.0e-6 500 1000
    fix 1 all nve # temp 300.0 300.0 1000  ### 
    fix 2 all langevin 285 285 1000 $(rand(1:100_000))) ### Tstart, Tstop, Dampening in fs  (1ps dampening in paper), seed
    fix 3 all print $(WriteOutFreq)  \"\$(step) \$(time) \$(temp) \$(press) \$(etotal) \$(pe) \$(ke) \$(evdwl) \$(ecoul) \$(epair) \$(ebond) \$(eangle) \$(edihed)\" file ./energy.dat screen no

    restart $(WriteOutFreq) ./sim.restart
    thermo $(WriteOutFreq)
    timestep 10.0
    run 1000000000")
    close(file)
end