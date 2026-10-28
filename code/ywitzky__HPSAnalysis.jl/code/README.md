# HPSAnalysis
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![CI](https://github.com/ywitzky/HPSAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ywitzky/HPSAnalysis.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/ywitzky/HPSAnalysis.jl/branch/main/graph/badge.svg?token=LHT7MZ2JP5)](https://codecov.io/gh/ywitzky/HPSAnalysis.jl)

This is a collection of scripts to set up molecular dynamics simulations using mostly the HOOMD simulation software for coarse-grained united residue models.

### Basic Usage 

```julia
###Initial coordinates of the proteins are calculated
(pos, Data) = HPS.CreateStartConfiguration(SimName,Path , Float32.([BoxLengthShort,BoxLengthShort*width_multiplier , BoxLengthShort]), Proteins, Sequences, Regenerate=false; Axis="y")#Erstellung der Start Conformation

###Inputs are written, including the simulation parameters, the interaction model, relevant dictionaries, and the start file
HPS.Setup.writeStartConfiguration("./$(protein)_slab","./$(SimName)_Start_slab.txt", Info, Sequences, BoxSize , 300_000_000, HOOMD=true, ; SimulationType="Calvados2" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH)


### simulate externally using "./src/Setup/Submit_HOOMD.py"

prot= HPS.initData("/Path/To/Data/"; Reparse=true, LoadAll=true, HOOMD=true)

prot.EquilibrationTime=...

### compute RG and it's auto correlation 
HPS.computeCOM!(prot)
HPS.computeRGSeries!(prot)
```


### Citation 
If you use HPSAnalysis, please cite the following papers: 

- Witzky, Yannick, Friederike Schmid, and Arash Nikoubashman. "From Heteropolymer Stiffness Distributions to Effective Homopolymers: A Conformational Analysis of Intrinsically Disordered Proteins." arXiv preprint arXiv:2504.11027 (2025).



