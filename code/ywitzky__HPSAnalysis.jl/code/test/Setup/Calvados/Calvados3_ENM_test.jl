using JSON

### test if all bonds are valid; not tested if all bonds that should be there are there....

@testset "Calvados3 ENM" begin 

BasePath = "$SetupTestPath/ENM_Test/"
if isdir(BasePath)
    rm(BasePath; force=true, recursive=true)
end
mkpath(BasePath)

DomainDict= Dict("RS31" => [[1,70], [90,150]], "RS31a" => [[1,75], [90,140]]) 
Proteins = ["RS31", "RS31", "RS31a"]
ProteinJSON= Dict("RS31" =>"$(PkgPath)/data/TestData/fold_rs31_full_data_0.json","RS31a" =>"$(PkgPath)/data/TestData/fold_rs31a_full_data_0.json" )
ProteinToCif= Dict("RS31" =>"$(PkgPath)/data/TestData/fold_rs31_model_0.cif","RS31a" =>"$(PkgPath)/data/TestData/fold_rs31a_model_0.cif" )


HPSAnalysis.RewriteCifToPDB(BasePath,ProteinToCif, Proteins )

for rcut in [7.0,8.0,10.0]
    ConstraintDict,Backbone_correction_Dict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = rcut, plDDTcut=90.0)
    RS31_r = getindex.(ConstraintDict["RS31"],3)
    RS31a_r = getindex.(ConstraintDict["RS31a"],3)
    @test all(RS31_r .<=rcut)
    @test all(RS31a_r .<=rcut) 
end


ConstraintDict,Backbone_correction_Dict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0, pae_cut=1.85)
RS31_i = getindex.(ConstraintDict["RS31"],1)
RS31_j = getindex.(ConstraintDict["RS31"],2)
dom1 = DomainDict["RS31"][1]
dom2 = DomainDict["RS31"][2]
 
isinbounds(x, dom1, dom2) = (x>=dom1[1]&& x<=dom1[2] )|| (x>=dom2[1]&& x<=dom2[2])

#check if ids are within bounds of folded regions
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31_i ))
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31_j ))


RS31a_i = getindex.(ConstraintDict["RS31a"],1)
RS31a_j = getindex.(ConstraintDict["RS31a"],2)
dom1 = DomainDict["RS31a"][1]
dom2 = DomainDict["RS31a"][2]

#check if ids are within bounds of folded regions
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31a_i ))
@test reduce(*, map(x -> isinbounds(x,dom1, dom2), RS31a_j ))


RS31_pae  =JSON.parsefile(ProteinJSON["RS31"])["pae"]
RS31a_pae =JSON.parsefile(ProteinJSON["RS31a"])["pae"]


### check if pae values are valid
check(i,j,_,_, pae, cut) =pae[i][j]<cut
@test reduce(*, map(x->check(x..., RS31_pae, 1.85),  ConstraintDict["RS31"]))
@test reduce(*, map(x->check(x..., RS31a_pae, 1.85),  ConstraintDict["RS31a"]))

CifData = Dict()
norm_RS31 = zeros(264)
RS31_plDDT = zeros(264)
RS31_x = zeros(264)
RS31_y = zeros(264)
RS31_z = zeros(264)
Masses = Dict("C"=>12.011, "N"=>14.007, "O"=>15.999, "S"=>32.065, "H"=>1.0078)
for line in split.(strip.(readlines(ProteinToCif["RS31"] )))
    if line[1]=="ATOM"
        cnt = parse(Int32, line[9])
        mass = Masses[line[3]]
        RS31_plDDT[cnt] += parse.(Float64,line[15])*mass
        RS31_x[cnt]+= parse.(Float64,line[11])*mass
        RS31_y[cnt]+= parse.(Float64,line[12])*mass
        RS31_z[cnt]+= parse.(Float64,line[13])*mass
        norm_RS31[cnt]  += mass
    end
end
RS31_plDDT ./= norm_RS31
RS31_x ./= norm_RS31
RS31_y ./= norm_RS31
RS31_z ./= norm_RS31


norm_RS31a = zeros(250)
RS31a_plDDT = zeros(250)
RS31a_x = zeros(250)
RS31a_y = zeros(250)
RS31a_z = zeros(250)
for line in split.(strip.(readlines(ProteinToCif["RS31a"] )))
    if line[1]=="ATOM"
        mass = Masses[line[3]]
        cnt = parse(Int32, line[9])
        RS31a_plDDT[cnt] += parse.(Float64,line[15])*mass
        RS31a_x[cnt]+= parse.(Float64,line[11])*mass
        RS31a_y[cnt]+= parse.(Float64,line[12])*mass
        RS31a_z[cnt]+= parse.(Float64,line[13])*mass
        norm_RS31a[cnt]  += mass
    end
end
RS31a_plDDT ./= norm_RS31a
RS31a_x ./= norm_RS31a
RS31a_y ./= norm_RS31a
RS31a_z ./= norm_RS31a

for (plDDTcut, pae_cut) in zip([80.0,90.0,92.0], [1.7, 1.85,2.0])
    local ConstraintDict,Backbone_correction_Dict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=plDDTcut, pae_cut=pae_cut)

    @test reduce(*, map(x->check(x..., RS31_pae ,pae_cut),  ConstraintDict["RS31"]))
    @test reduce(*, map(x->check(x..., RS31a_pae,pae_cut),  ConstraintDict["RS31a"]))

    local RS31_i  = getindex.(ConstraintDict["RS31"] ,1)
    local RS31_j  = getindex.(ConstraintDict["RS31"] ,2)
    local RS31a_i = getindex.(ConstraintDict["RS31a"],1)
    local RS31a_j = getindex.(ConstraintDict["RS31a"],2)
    
    @test reduce(*, map(i-> RS31_plDDT[i] >plDDTcut, RS31_i ))
    @test reduce(*, map(i-> RS31_plDDT[i] >plDDTcut, RS31_j ))
    @test reduce(*, map(i-> RS31a_plDDT[i]>plDDTcut, RS31a_i))
    @test reduce(*, map(i-> RS31a_plDDT[i]>plDDTcut, RS31a_j))
end


###test whether all bonds are found 
DomainDict= Dict("RS31" => [[1,70], [90,130]], "RS31a" => [[1,60], [90,140]]) 
ConstraintDict,Backbone_correction_Dict = HPSAnalysis.Setup.DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA",precision=10^12)


valid_plDDT = [i for (i, val) in  enumerate(RS31_plDDT[1:70]) if val>90 ]
RS31_dom1_pairs = [(i,j) for i in valid_plDDT for j in valid_plDDT if i<j-2 && RS31_pae[i][j]<1.85 && RS31_pae[j][i]<1.85 ]

valid_plDDT = [i+89 for (i, val) in  enumerate(RS31_plDDT[90:130]) if val>90 ]
RS31_dom2_pairs = [(i,j) for i in valid_plDDT for j in valid_plDDT if i<j-2 && RS31_pae[i][j]<1.85  && RS31_pae[j][i]<1.85]

valid_plDDT = [i for (i, val) in  enumerate(RS31a_plDDT[1:60]) if val>90 ] 
RS31a_dom1_pairs = [(i,j) for i in valid_plDDT for j in valid_plDDT if i<j-2 && RS31a_pae[i][j]<1.85 && RS31a_pae[j][i]<1.85]

valid_plDDT = [i+89 for (i, val) in  enumerate(RS31a_plDDT[90:140]) if val>90 ]
RS31a_dom2_pairs = [(i,j) for i in valid_plDDT for j in valid_plDDT if i<j-2 && RS31a_pae[i][j]<1.85  && RS31a_pae[j][i]<1.85]


RS31_dist  = [ (RS31_x[i]  -RS31_x[j])^2+(RS31_y[i]  -RS31_y[j])^2+(RS31_z[i]  -RS31_z[j])^2 for i in 1:150, j in 1:150 ]
RS31a_dist = [ (RS31a_x[i]-RS31a_x[j])^2+(RS31a_y[i]-RS31a_y[j])^2+(RS31a_z[i]-RS31a_z[j])^2 for i in 1:150, j in 1:150 ]
RS31_dist  .= sqrt.(RS31_dist)
RS31a_dist .= sqrt.(RS31a_dist)

RS31_pairs  = [(i,j,RS31_dist[i,j]/10.0)  for (i,j) in vcat(RS31_dom1_pairs ,RS31_dom2_pairs ) if RS31_dist[i,j]<9.0 ]
RS31a_pairs = [(i,j,RS31a_dist[i,j]/10.0) for (i,j) in vcat(RS31a_dom1_pairs,RS31a_dom2_pairs) if RS31a_dist[i,j]<9.0 ]

@test all(map((x,y) -> x[1]==y[1]&& x[2]==y[2] && x[3]≈y[3], RS31_pairs ,  ConstraintDict["RS31"]))
@test all(map((x,y) -> x[1]==y[1]&& x[2]==y[2] && x[3]≈y[3], RS31a_pairs,  ConstraintDict["RS31a"]))

### test the generation of unfolded regions
Proteins = ["RS31", "RS31", "RS31a"]
Sequences = [HPSAnalysis.ProteinSequences.NameToSeq[prot] for prot in Proteins]
UnfoldedRegions = HPSAnalysis.Setup.GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)

@test UnfoldedRegions == Dict("RS31a" => [(60, 90), (140, 250)], "RS31" => [(70, 90), (130, 264)])


### check if all backbones are fully existent by automatic build

Sim = HPSAnalysis.SimData()
Sim.BasePath = SetupTestPath
mkpath("$(Sim.BasePath)/InitFiles/CifFiles/")
for Protein in Set(Proteins)
    cp(ProteinToCif[Protein], "$(Sim.BasePath)/InitFiles/CifFiles/$(Protein).cif"; force=true)
end
(NBonds, B_types, B_typeid, B_groups, harmonic ) = HPSAnalysis.Setup.BuildENMModel(Sim, DomainDict, Proteins, Sequences, ProteinJSON)


offsets = vcat([0], cumsum(length.(Sequences))[1:end-1])
indices =  vcat([collect(1+off:length(seq)-1+off) for (seq, off) in zip(Sequences, offsets)]...) .- 1 ### use c indexing

@test all(map(i-> (i,i+1) in B_groups, indices )) ## -1 because of python/c++ indexing

### test read and write procedure
HPSAnalysis.Setup.WriteENM_HOOMD_Indices("$SetupTestPath/HOOMD_Setup/ENM_indices.txt", (NBonds, B_types, B_typeid, B_groups, harmonic))

data  = [strip.(split(line, ",")) for line in readlines("$SetupTestPath/HOOMD_Setup/ENM_indices.txt")[2:end]]
ids = parse.(Int32, getindex.(data,3))
types  = getindex.(data,2)

@test [val != "O-O" ? parse(Int32,split(val,"_")[2]) : 0 for val in types] == ids ### test if id matches the one in type name
@test ["$(harmonic[v])" for v in types] == ["$n, $m" for (n,m) in zip(getindex.(data, 6), getindex.(data, 7))]
if PythonTests
    (NBonds_read, B_types_read, B_typeid_read, B_groups_read, harmonic_read) = sim.read_ENM_HOOD_indices("$SetupTestPath/HOOMD_Setup/ENM_indices.txt")

    @test NBonds  == NBonds_read
    @test all(B_types  == B_types_read)
    @test B_typeid == B_typeid_read
    @test B_groups == B_groups_read
    @test all(map(key-> (harmonic[key][:k]== harmonic_read[key]["k"])&& (harmonic[key][:r]≈harmonic_read[key]["r"]), collect(keys(harmonic))))
end


### Test using Constraints for copies of the same protein
### Backbone_correction_Dict should contain all backbones connection in folded regions but they dont in the test
ConstraintDict = Dict("Prot1"=> [(1,4,0.5,3), (10,12, 1.2,4), (14,15,0.2,5)], "Prot2"=> [(3,7,0.5,6), (3,9, 1.3,7), (6,9,0.33,8)])
Backbone_correction_Dict = Dict("Prot1"=> [(2,3,0.37,1), (4,5,0.38,0)], "Prot2"=> [(5,6,0.5,2), (6,7,0.38,0)])
DomainDict = Dict("Prot1" => [(1,4), (10,15)], "Prot2"=> [(3,9)])
Seq_Dict = Dict("Prot1"=> "ABCABCABCABCABC", "Prot2"=> "ABABABABAB")
Proteins = ["Prot1", "Prot1", "Prot2", "Prot1", "Prot2"]
Sequences = [Seq_Dict[x] for x in Proteins]



### check generation of bonds for IDRs
UnfoldedRegions = HPSAnalysis.Setup.GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)
(NBonds, B_types, B_typeid, B_groups, _) = HPSAnalysis.Setup.CombineBackboneAndENM(Proteins, Sequences, (0,["O-O"], [], [], Dict()), UnfoldedRegions, Backbone_correction_Dict)

@test NBonds == 21
@test B_types == ["O-O"]
@test B_typeid == zeros(21)
UnfoldedGroups = sort(vcat([[(4,5).+x , (5,6).+x , (6,7).+x ,(7,8).+x , (8,9).+x] for x in [0, 15, 40]]..., [[(0,1).+x , (1,2).+x , (8,9).+x ] for x in [30, 55]]...))  ### use c indexing

@test B_groups == UnfoldedGroups

### only generate additional bonds due to enm
(NBonds, B_types, B_typeid, B_groups, harmonic) = HPSAnalysis.Setup.ComputeHOOMD_ENMIndices(ConstraintDict, Backbone_correction_Dict, Sequences, Proteins)

FoldedGroups= [(1, 2), (3, 4), (16, 17), (18, 19), (34, 35), (35, 36), (41, 42), (43, 44), (59, 60), (60, 61), (0, 3), (9, 11), (13, 14), (15, 18), (24, 26), (28, 29), (32, 36), (32, 38), (35, 38), (40, 43), (49, 51), (53, 54), (57, 61), (57, 63), (60, 63)]

@test NBonds == 25
@test all(B_types .== ["O-O", "BB_1", "BB_2", "ENM_3", "ENM_4", "ENM_5", "ENM_6","ENM_7", "ENM_8"])
@test all(B_typeid .== [1,0,1,0,2,0,1,0,2,0, 3,4,5,3,4,5,6,7,8,3,4,5,6,7,8]) # has a shift because default BB is 1
@test all(B_groups .== FoldedGroups) # shift from 1 -> 0, because of julia -> python


harmonic_test = Dict{String, Dict{Symbol, Float64}}()
cnt = 5
for Protein in Set(Proteins)
    for (i,j,r0, ind) = ConstraintDict[Protein]
        bondname = "ENM_$(ind)"
        harmonic_test[bondname] = Dict(:r => r0, :k => 700)
        cnt+= 1
    end
end

harmonic_test["BB_1"] = Dict(:k => 8033.0, :r => 0.37)
harmonic_test["BB_2"] = Dict(:k => 8033.0, :r => 0.5 )
harmonic_test["O-O"]  = Dict(:k => 8033.0, :r => 0.38 )

@test harmonic_test == harmonic


### combine ENM and backbone of IDRs 

(Comb_Bonds, Comb_types, Comb_typeid, Comb_groups, Comb_harmonic) = HPSAnalysis.Setup.CombineBackboneAndENM(Proteins, Sequences, (NBonds, B_types, B_typeid, B_groups, harmonic), UnfoldedRegions, Backbone_correction_Dict)

@test Comb_Bonds == 25+21
@test Comb_types == ["O-O", "BB_1", "BB_2", "ENM_3", "ENM_4", "ENM_5", "ENM_6","ENM_7", "ENM_8"]
@test Comb_typeid == Int32.(vcat(zeros(21),  [1,0,1,0,2,0,1,0,2,0, 3,4,5,3,4,5,6,7,8,3,4,5,6,7,8]))
@test Comb_groups == vcat(UnfoldedGroups, FoldedGroups)



### test write out and reading again
filepath = BasePath*"/ENM_indices.txt"
HPSAnalysis.Setup.WriteENM_HOOMD_Indices(filepath, (Comb_Bonds, Comb_types, Comb_typeid, Comb_groups, Comb_harmonic))

data  = [strip.(split(line, ",")) for line in readlines(filepath)[2:end]]
ids = parse.(Int32, getindex.(data,3))
types  = getindex.(data,2)

@test [val != "O-O" ? parse(Int32,split(val,"_")[2]) : 0 for val in types] == ids ### test if id matches the one in type name
@test ["$(Comb_harmonic[v])" for v in types] == ["$n, $m" for (n,m) in zip(getindex.(data, 6), getindex.(data, 7))]

testpos = [ [1 0 0;], [0 2 0;], [0 0 3], [1 2 3;]]

### rotate 90 degrees around x, rotate 180 around x and 270 around y
RotMats= [ [1 0 0; 0 0 1; 0 -1 0],  [0 0 -1; 0 -1 0; -1 0 0]]

solution = [[1 0 0;], [0 0 -1;], [0 0 2;], [0 -2 0;], [0 -3 0;], [-3 0 0;], [1 -3 2;], [-3 -2 -1]]

test = [pos*mat for pos in testpos for mat in RotMats[1:1]]
test_vals= [HPSAnalysis.Setup.rotate_randomly(pos, RotMats[1:1]) for pos in testpos]
@test all(test_vals.≈solution[1:2:end])

test = [pos*mat for pos in testpos for mat in RotMats[2:2]]
test_vals= [HPSAnalysis.Setup.rotate_randomly(pos, RotMats[2:2]) for pos in testpos]
@test all(test_vals.≈solution[2:2:end])
end