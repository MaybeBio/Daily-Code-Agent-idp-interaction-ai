@testset "MeanPairEnergy" begin

Sim = HPSAnalysis.SimData()

Sim.BasePath= TestPath
mkpath("$(TestPath)/HOOMD_Setup/")
### test parameter reader before manipulating it for test purposes; data in files is always in nm wheres analysis works in Angstroem -> conversion
function write_param(yu_pref, κ,rcut_yu, rcut_ah)
    file=open("$(TestPath)/HOOMD_Setup/Params.txt","w")
    write(file,"yk_prefactor: $(yu_pref/10.0)\n")
    write(file,"kappa: $(κ*10.0)\n")
    write(file,"YukawaCutoff: $(rcut_yu/10.0)\n")
    write(file,"AHCutoff: $(rcut_ah/10.0)\n")
    close(file)
end

test = []
for pref in 0.5:1.5
    for κ_test in 10.0:2.0:12.0
        for r_yu in 35:5:40
            for r_ah in 18:2:20
                write_param(pref, κ_test, r_yu,r_ah)

                yu_pref, κ, rcut_yu, rcut_ah = HPSAnalysis.read_FF_Interaction_Parameters(Sim)

                push!(test, yu_pref≈ pref && κ_test≈κ && r_yu≈rcut_yu && r_ah ≈ rcut_ah)
            end
        end
    end
end

@test all(test) ### report all reads as one test


### Actual test of the functions

# Build a tiny test system where PBC matters
#Random.seed!(42)                       # reproducibility

L = 100.0                                 # orthorhombic box side length
Sim = HPSAnalysis.SimData()
Sim.NAtoms      = 6
Sim.NChains     = 2
Sim.NSteps      = 1                     # **single frame**
Sim.ChainStart  = [1,4]
Sim.ChainStop   = [3,6]
Sim.ChainLength = [3,3]
Sim.x           = zeros(6,1)
Sim.y           = zeros(6,1)
Sim.z           = zeros(6,1)
Sim.BoxLength   = [L, L, L]
Sim.Charges     = zeros(6)
Sim.IDs         = [1,2,3,1,2,3]
Sim.IDToResName = Dict(1=>"A",2=>"B",3=>"C")
Sim.IDToSigmas  = collect(5:7)
Sim.IDToLambdas = collect(LinRange(0.2,1,3))
Sim.Masses      = ones(Sim.NAtoms)
Sim.ChainMasses = [3,3]
Sim.EquilibrationTime = 1 
Sim.RGMeasureStep     = 1

NTypes= length(keys(Sim.IDToResName))
# Coordinates – place one pair across the box boundary
# Chain 1
Sim.x[1,1] = 0.10;  Sim.y[1,1] = 0.00;  Sim.z[1,1] = 0.00   # residue 1
Sim.x[2,1] = 0.90;  Sim.y[2,1] = 0.00;  Sim.z[2,1] = 0.00   # residue 2 (near right wall)
Sim.x[3,1] = 0.50;  Sim.y[3,1] = 0.50;  Sim.z[3,1] = 0.00   # residue 3

# Chain 2
Sim.x[4,1] = 0.20;  Sim.y[4,1] = 0.00;  Sim.z[4,1] = 0.00   # residue 4 (will be close to #2 across PBC)
Sim.x[5,1] = 0.80;  Sim.y[5,1] = 0.50;  Sim.z[5,1] = 0.00   # residue 5
Sim.x[6,1] = 0.40;  Sim.y[6,1] = 0.50;  Sim.z[6,1] = 0.00   # residue 6

Sim.x .*= L
Sim.y .*= L
Sim.z .*= L

Sim.x_uw= Sim.x
Sim.y_uw= Sim.y
Sim.z_uw= Sim.z


# Charges – only a few residues are charged 
Sim.Charges[1] =  +1.0     # residue 1
Sim.Charges[2] =  -1.0     # residue 2 (direct neighbour of #1 → should be *ignored*)
Sim.Charges[4] =  +1.0     # residue 4 (will interact with #2 across PBC)
Sim.Charges[5] =  -1.0     # residue 4 (will interact with #2 across PBC)


# Manual reference calculation (exactly the same maths as the routine)
kappa = 1.0/2.222
yu_pref = 4.61008
ϵ_ah = 0.89
rcut_yu = 35.0
rcut_ah = 15.0
yu(r,q1,q2) = yu_pref * q1 * q2 / r * exp(-kappa * r)
ah(r,σ,λ) = HPSAnalysis.ashbaugh_hatch_potential(Float64(r), σ,λ, ϵ_ah)

### manipulate parameter reader to return test parameter
function HPSAnalysis.read_FF_Interaction_Parameters(Sim::HPSAnalysis.SimData{R,I}) where {R<:Real, I<:Integer}
    return R(yu_pref), R(kappa), R(rcut_yu), R(rcut_ah)
end

@test all([isapprox(HPSAnalysis.ashbaugh_hatch_potential_fast(r,σ,λ,ϵ),HPSAnalysis.ashbaugh_hatch_potential(r,σ,λ,ϵ),atol=10^-12,rtol=10^-12) for σ in 2:0.3:5.0 for r in 0.8*σ:0.1:2.0*σ for λ in 0:0.2:1 for ϵ in 0.0:0.5:5.0])

# minimum‑image helper
min_image(d, L) = abs(d) > L/2 ? d - sign(d)*L : d

function pair_energy(i,j)
    dx = Sim.x[i,1] - Sim.x[j,1]
    dy = Sim.y[i,1] - Sim.y[j,1]
    dz = Sim.z[i,1] - Sim.z[j,1]
    dx = min_image(dx, L)
    dy = min_image(dy, L)
    dz = min_image(dz, L)
    r  = sqrt(dx*dx + dy*dy + dz*dz)
    σ = (Sim.IDToSigmas[ Sim.IDs[i]]+Sim.IDToSigmas[ Sim.IDs[j]])/2.0
    λ = (Sim.IDToLambdas[Sim.IDs[i]]+Sim.IDToLambdas[Sim.IDs[j]])/2.0

    en_yu = r < rcut_yu ? yu(r, Sim.Charges[i], Sim.Charges[j]) : 0.0
    en_ah = r < rcut_ah ? ah(r,σ,λ) : 0.0
    return en_yu , en_ah
end

chainLen = Sim.ChainLength[1]

# reference matrices (the same shape as the fast routine)
yu_inter   = zeros(Float64, chainLen, chainLen)
yu_interAA = zeros(Float64, NTypes, NTypes)

ah_inter   = zeros(Float64, chainLen, chainLen)
ah_interAA = zeros(Float64, NTypes, NTypes)

# loop over chain pairs
for C1 in 1:Sim.NChains
    for C2 in C1+1:Sim.NChains
        for (r1,i) in enumerate(Sim.ChainStart[C1]:Sim.ChainStop[C1])
            for (r2,j) in enumerate(Sim.ChainStart[C2]:Sim.ChainStop[C2])
            # ---- compute the energy for this pair ----
            en_yu, en_ah = pair_energy(i,j)

            yu_inter[r1,r2]   += en_yu
            ah_inter[r1,r2]   += en_ah

            yu_interAA[Sim.IDs[i],Sim.IDs[j]] += en_yu
            ah_interAA[Sim.IDs[i],Sim.IDs[j]] += en_ah
            end
        end
    end
end

# Symmetrise and halve the diagonal (exactly what the routine does)
function symmetrise!(M)
    M .+= transpose(M)
    for i in 1:size(M,1)
        M[i,i] *= 0.5
    end
    M
end

symmetrise!(yu_inter)
symmetrise!(ah_inter)

symmetrise!(yu_interAA)
symmetrise!(ah_interAA)

# Call the fast routine for single frame, bounds large enough to not exclude anything
HPSAnalysis.computeCOM!(Sim)    ### needs to be precomputed
HPSAnalysis.computeClustersByChainCOM(Sim; Cutoff=70.0) ### needs to be precomputed
bounds = Vector{Tuple{Float32, Float32}}([Float32.((-L/2.0, L/2.0)), Float32.((-L/2.0, L/2.0)),Float32.((-L/2.0, L/2.0))])
yu_inter_fast, yu_interAA_fast, ah_inter_fast, ah_interAA_fast = HPSAnalysis.computeMeanPairEnergyMatrix_naiv(Sim, bounds; ϵ_ah=Float32(ϵ_ah))


AAcount_per_chain =[1,1,1]
pair_matrix = [i == j  ?  
        AAcount_per_chain[i] * AAcount_per_chain[i]*  0.5 : 
        AAcount_per_chain[i] * AAcount_per_chain[j] 
        for i in 1:NTypes, j in 1:NTypes
    ]
    
yu_interAA ./= pair_matrix
ah_interAA ./= pair_matrix


# Compare the two results
@test all(isapprox.(yu_inter, yu_inter_fast; rtol=1e-6, atol=1e-6))
@test all(isapprox.(ah_inter, ah_inter_fast; rtol=1e-6, atol=1e-6))

@test all(isapprox.(yu_interAA, yu_interAA_fast; rtol=1e-6, atol=1e-4))
@test all(isapprox.(ah_interAA, ah_interAA_fast; rtol=1e-6, atol=1e-4))




######### Now more steps and bounds

L = 100.0
Sim = HPSAnalysis.SimData()
Sim.NAtoms      = 12
Sim.NChains     = 3
Sim.NSteps      = 3                     # Multiple time steps
Sim.ChainStart  = [1, 5, 9]
Sim.ChainStop   = [4, 8, 12]
Sim.ChainLength = [4, 4, 4]
Sim.x           = zeros(12, 3)          # 12 atoms × 3 steps
Sim.y           = zeros(12, 3)
Sim.z           = zeros(12, 3)
Sim.BoxLength   = [L, L, L]
Sim.Charges     = zeros(12)
Sim.IDs         = [1,2,3,3,1,2,3,3,1,2,3,3]
Sim.IDToResName = Dict(i => string('A' + i - 1) for i in 1:maximum(Sim.IDs))  # Unique residue names
Sim.IDToSigmas  = [3,4,5]
Sim.IDToLambdas = collect(LinRange(0.2, 1, 3))
Sim.Masses      = ones(12)
Sim.ChainMasses = [4, 4, 4]
Sim.EquilibrationTime = 1
Sim.RGMeasureStep     = 1

IDToCharges = [1.0,0,-1.0]
Sim.Charges = [IDToCharges[id] for id in Sim.IDs]

# Coordinates setup (all in [0,1] range before scaling)
# Step 1: All chains fully within bounds (x ∈ [20,80])
# Step 2: Chain 1 partially outside (x=10 < 20)
# Step 3: Chains 1 & 2 partially outside (x=10 and x=90)
for step in 1:3
    # Chain 1 (atoms 1-4)
    Sim.x[1:4, step] = [30, 40, 50, 60] 
    Sim.y[1:4, step] .= 50 
    Sim.z[1:4, step] .= 47.5
    
    # Chain 2 (atoms 5-8)
    Sim.x[5:8, step] = [30, 40, 50, 60] 
    Sim.y[5:8, step] .= 55
    Sim.z[5:8, step] .= 50 
    
    # Chain 3 (atoms 9-12)
    Sim.x[9:12, step] = [30, 40, 50, 60] 
    Sim.y[9:12, step] .= 45
    Sim.z[9:12, step] .= 52.5 
end

# Modify specific positions to create boundary conditions
Sim.x[1, 2] = 5    # Chain1 atom1 outside bounds in step2 (x=10 < 20)
Sim.x[5, 3] = 95   # Chain2 atom5 outside bounds in step3 (x=90 > 80)

# center box
Sim.x .-= L/2
Sim.y .-= L/2
Sim.z .-= L/2

Sim.x_uw = copy(Sim.x)
Sim.y_uw = copy(Sim.y)
Sim.z_uw = copy(Sim.z)

# Reference calculation setup
NTypes = length(keys(Sim.IDToResName))
chainLen = maximum(Sim.ChainLength)

# Initialize accumulators for mean calculation
yu_inter_acc   = zeros(Float64, chainLen, chainLen)
yu_interAA_acc = zeros(Float64, NTypes, NTypes)
ah_inter_acc   = zeros(Float64, chainLen, chainLen)
ah_interAA_acc = zeros(Float64, NTypes, NTypes)

# Define bounds: only consider chains with ALL atoms in x∈[-40,40], y/z∈[-50,50]
bounds = [(-40.0f0, 40.0f0), (-50.0f0, 50.0f0),  (-50.0f0, 50.0f0)]

Sim.EquilibrationTime=1
Sim.NSteps=1
NValid=0
# Manual reference calculation over valid steps
for step in Sim.EquilibrationTime:Sim.NSteps
    # Check which chains are fully within bounds
    valid_chains = Int[]
    for chain in 1:Sim.NChains
        start_idx = Sim.ChainStart[chain]
        stop_idx  = Sim.ChainStop[chain]
        valid = true
        
        for atom in start_idx:stop_idx
            # Check all 3 dimensions
            if !(bounds[1][1] <= Sim.x[atom, step] <= bounds[1][2]) ||
               !(bounds[2][1] <= Sim.y[atom, step] <= bounds[2][2]) ||
               !(bounds[3][1] <= Sim.z[atom, step] <= bounds[3][2])
                valid = false
                break
            end
        end
        
        if valid
            push!(valid_chains, chain)
        end
    end
    
    # Skip step if fewer than 2 valid chains
    #length(valid_chains) < 2 && continue
    NValid += length(valid_chains)
    
    # Process valid chain pairs
    for chain1 in valid_chains
        for chain2 in 1:Sim.NChains
            if chain1==chain2 continue end
            
            for (r1, i) in enumerate(Sim.ChainStart[chain1]:Sim.ChainStop[chain1])
                for (r2, j) in enumerate(Sim.ChainStart[chain2]:Sim.ChainStop[chain2])
                    # Compute pair distance with PBC
                    dx = min_image(Sim.x[i,step] - Sim.x[j,step], L)
                    dy = min_image(Sim.y[i,step] - Sim.y[j,step], L)
                    dz = min_image(Sim.z[i,step] - Sim.z[j,step], L)
                    r = sqrt(dx^2 + dy^2 + dz^2)
                    
                    # Get interaction parameters
                    id_i = Sim.IDs[i]
                    id_j = Sim.IDs[j]
                    σ = (Sim.IDToSigmas[ id_i] + Sim.IDToSigmas[ id_j]) / 2.0
                    λ = (Sim.IDToLambdas[id_i] + Sim.IDToLambdas[id_j]) / 2.0
                    
                    # Compute energies
                    en_yu = r < rcut_yu ? yu(r,Sim.Charges[i],Sim.Charges[j]) : 0.0
                    en_ah = r < rcut_ah ? ah(r, σ, λ) : 0.0

                    # Accumulate chain-chain interactions
                    yu_inter_acc[r1, r2] += en_yu
                    ah_inter_acc[r1, r2] += en_ah
                    
                    # Accumulate residue-type interactions
                    type_i = Sim.IDs[i]
                    type_j = Sim.IDs[j]

                    yu_interAA_acc[type_i, type_j] += en_yu
                    ah_interAA_acc[type_i, type_j] += en_ah
                end
            end
        end
    end
end

# Compute means (only over valid steps)
if NValid > 0
    yu_inter_acc   ./= NValid
    ah_inter_acc   ./= NValid
    yu_interAA_acc ./= NValid
    ah_interAA_acc ./= NValid
end

symmetrise!(yu_inter_acc)
symmetrise!(ah_inter_acc)
symmetrise!(yu_interAA_acc)
symmetrise!(ah_interAA_acc)

# Precompute required properties
HPSAnalysis.computeCOM!(Sim)
HPSAnalysis.computeClustersByChainCOM(Sim; Cutoff=200.0)

# Compute via target function (with bounds filtering)
yu_inter_fast, yu_interAA_fast, ah_inter_fast, ah_interAA_fast = 
    HPSAnalysis.computeMeanPairEnergyMatrix_naiv(Sim, bounds; ϵ_ah=Float32(ϵ_ah))


AAcount_per_chain =[1,1,2]
pair_matrix = [i == j  ?  
        AAcount_per_chain[i] * AAcount_per_chain[i] *  0.5 : 
        AAcount_per_chain[i] * AAcount_per_chain[j] 
        for i in 1:NTypes, j in 1:NTypes
    ]

ah_interAA_acc ./= pair_matrix 
yu_interAA_acc   ./= pair_matrix

@test all(isapprox.(yu_inter_acc, yu_inter_fast; rtol=1e-6, atol=1e-6))
@test all(isapprox.(ah_inter_acc, ah_inter_fast; rtol=1e-6, atol=1e-6))
@test all(isapprox.(yu_interAA_acc, yu_interAA_fast; rtol=1e-6, atol=1e-2))
@test all(isapprox.(ah_interAA_acc, ah_interAA_fast; rtol=1e-6, atol=1e-2))


end