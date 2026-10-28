
# -------------------------------------------------------------------------
#   Helper: mimic the Yukawa and Ashbaugh‑Hatch potentials used by the package
# -------------------------------------------------------------------------
yu(r, q1, q2) = yu_pref * q1 * q2 / r * exp(-kappa * r)
ah(r, σ, λ) = HPSAnalysis.ashbaugh_hatch_potential(Float64(r), σ, λ, ϵ_ah)

# -------------------------------------------------------------------------
#   Minimal helper: periodic‑image distance (same implementation as the
#   reference code in your original test)
# -------------------------------------------------------------------------
min_image(d, L) = abs(d) > L/2 ? d - sign(d) * L : d

# -------------------------------------------------------------------------
#   Fake interaction parameters – we force the routine to use these values
#   by monkey‑patching `read_FF_Interaction_Parameters`.
# -------------------------------------------------------------------------
const yu_pref = 4.61008
const kappa   = 1.0/2.222
const rcut_yu = 35.0
const rcut_ah = 15.0
const ϵ_ah    = 0.89

function HPSAnalysis.read_FF_Interaction_Parameters(::HPSAnalysis.SimData{R,I}) where {R<:Real,I<:Integer}
    return R(yu_pref), R(kappa), R(rcut_yu), R(rcut_ah)
end

# -------------------------------------------------------------------------
#   Monkey‑patch `getHOOMDBody` so that it returns a user‑defined vector.
#   (The original implementation reads a GSD file – we want to stay completely
#   in‑memory for the test.)
# -------------------------------------------------------------------------
#   The body vector is stored in a global constant that the test can modify.
const _test_body_vec = Ref{Vector{Int}}()

import HPSAnalysis: getHOOMDBody
function getHOOMDBody(::HPSAnalysis.SimData{R,I}) where {R<:Real,I<:Integer}
    return _test_body_vec[]
end

# symmetrise and halve the diagonal (exactly what the library routine does)
function symmetrise!(M)
    M .+= transpose(M)
    for i in 1:size(M,1)
        M[i,i] *= 0.5
    end
    M
end

# =============================================================================
#   1️⃣  SINGLE‑STEP INTRA‑CHAIN TEST
# =============================================================================
@testset "IntraChainEnergy – single step" begin
    # --------------------------------------------------------------
    #   Build a tiny simulation with ONE chain (6 residues)
    # --------------------------------------------------------------
    L = 100.0                                   # orthorhombic box side length
    Sim = HPSAnalysis.SimData()
    Sim.NAtoms      = 6
    Sim.NChains     = 1
    Sim.NSteps      = 1                         # only one frame
    Sim.ChainStart  = [1]
    Sim.ChainStop   = [6]
    Sim.ChainLength = [6]
    Sim.x           = zeros(6, 1)                # positions will be overwritten below
    Sim.y           = zeros(6, 1)
    Sim.z           = zeros(6, 1)
    Sim.BoxLength   = [L, L, L]
    Sim.Charges     = zeros(6)
    Sim.IDs         = [1,2,3,1,2,3]            # three residue types (A,B,C)
    Sim.IDToResName = Dict(1=>"A", 2=>"B", 3=>"C")
    Sim.IDToSigmas  = collect(5:7)               # σA=5, σB=6, σC=7
    Sim.IDToLambdas = collect(LinRange(0.2,1,3)) # λA=0.2, λB=0.6, λC=1.0
    Sim.Masses      = ones(6)
    Sim.ChainMasses = [6]
    Sim.EquilibrationTime = 1
    Sim.RGMeasureStep     = 1

    # -----------------------------------------------------------------
    #   Place residues so that a few pairs cross the periodic boundary
    # -----------------------------------------------------------------
    Sim.x[1,1] = 0.10;  Sim.y[1,1] = 0.00;  Sim.z[1,1] = 0.00   # residue 1 (type A)
    Sim.x[2,1] = 0.90;  Sim.y[2,1] = 0.00;  Sim.z[2,1] = 0.00   # residue 2 (type B) – near right wall
    Sim.x[3,1] = 0.55;  Sim.y[3,1] = 0.50;  Sim.z[3,1] = 0.00   # residue 3 (type C)
    Sim.x[4,1] = 0.30;  Sim.y[4,1] = 0.20;  Sim.z[4,1] = 0.00   # residue 4 (type A)
    Sim.x[5,1] = 0.80;  Sim.y[5,1] = 0.70;  Sim.z[5,1] = 0.00   # residue 5 (type B)
    Sim.x[6,1] = 0.40;  Sim.y[6,1] = 0.55;  Sim.z[6,1] = 0.00   # residue 6 (type C)

    # Convert from reduced coordinates [0,1] to Å
    Sim.x .*= L; Sim.y .*= L; Sim.z .*= L

    # copy to unwrapped coordinates
    Sim.x_uw= Sim.x
    Sim.y_uw= Sim.y
    Sim.z_uw= Sim.z

    # -----------------------------------------------------------------
    #   Charge pattern – we give a non‑zero charge to a few residues
    # -----------------------------------------------------------------
    Sim.Charges[1] =  +1.0   # A
    Sim.Charges[2] =  -1.0   # B (direct neighbour → should be ignored if beyond cutoff)
    Sim.Charges[5] =  +1.0   # B (will interact with residue 2 across the box)

    # -----------------------------------------------------------------
    #   Body vector: particles with body == -1 are *active*.
    #   Here residues 2,4,5,6 are active, residues 1 and 3 are inert.
    # -----------------------------------------------------------------
    _test_body_vec[] = [0, -1, 0, -1, -1, -1]

    # -----------------------------------------------------------------
    #   Reference calculation (hand‑written)
    # -----------------------------------------------------------------
    chain_len = Sim.ChainLength[1]
    NTypes    = length(keys(Sim.IDToResName))

    yu_ref   = zeros(Float64, chain_len, chain_len)
    ah_ref   = zeros(Float64, chain_len, chain_len)
    yu_refAA = zeros(Float64, NTypes, NTypes)
    ah_refAA = zeros(Float64, NTypes, NTypes)

    # loop over *all* unordered pairs inside the single chain
    for (i_rel, i) in enumerate(Sim.ChainStart[1]:Sim.ChainStop[1])
        for (j_rel, j) in enumerate(Sim.ChainStart[1]:Sim.ChainStop[1])
            i_rel >= j_rel && continue                     # skip self‑interaction
            # ---- body mask -------------------------------------------------
            !(( _test_body_vec[][i] == -1 ) && ( _test_body_vec[][j] == -1 )) && continue
            # ---- minimum image distance ------------------------------------
            dx = min_image(Sim.x[i,1] - Sim.x[j,1], L)
            dy = min_image(Sim.y[i,1] - Sim.y[j,1], L)
            dz = min_image(Sim.z[i,1] - Sim.z[j,1], L)
            r  = sqrt(dx*dx + dy*dy + dz*dz)

            # ---- interaction parameters ------------------------------------
            σ = (Sim.IDToSigmas[ Sim.IDs[i] ] + Sim.IDToSigmas[ Sim.IDs[j] ]) / 2.0
            λ = (Sim.IDToLambdas[ Sim.IDs[i] ] + Sim.IDToLambdas[ Sim.IDs[j] ]) / 2.0

            # ---- Yukawa ----------------------------------------------------
            en_yu = r < rcut_yu ? yu(r, Sim.Charges[i], Sim.Charges[j]) : 0.0
            # ---- Ashbaugh‑Hatch --------------------------------------------
            en_ah = r < rcut_ah ? ah(r, σ, λ) : 0.0

            # ---- accumulate ------------------------------------------------
            yu_ref[i_rel, j_rel]   += en_yu
            ah_ref[i_rel, j_rel]   += en_ah
            yu_refAA[ Sim.IDs[i], Sim.IDs[j] ] += en_yu
            ah_refAA[ Sim.IDs[i], Sim.IDs[j] ] += en_ah
        end
    end

    symmetrise!(yu_ref);   symmetrise!(ah_ref)
    symmetrise!(yu_refAA); symmetrise!(ah_refAA)

    # -----------------------------------------------------------------
    #   Call the routine that should *only* see residues with body == -1
    # -----------------------------------------------------------------
    bounds = Vector{Tuple{Float32,Float32}}(undef,3)
    bounds[1] = (Float32(-L/2.0), Float32(L/2.0))
    bounds[2] = (Float32(-L/2.0), Float32(L/2.0))
    bounds[3] = (Float32(-L/2.0), Float32(L/2.0))

    HPSAnalysis.computeCOM!(Sim)    ### needs to be precomputed
    HPSAnalysis.computeClustersByChainCOM(Sim; Cutoff=70.0) ### needs to be precomputed
    # NOTE: `getIntraChainEnergyMap` is assumed to return a 4‑tuple of matrices
    # exactly like `computeMeanPairEnergyMatrix_naiv` does.
    yu_fast, yuAA_fast, ah_fast, ahAA_fast = HPSAnalysis.getIntraChainEnergyMap(Sim, bounds; ϵ_ah=Float32(ϵ_ah))

    # -----------------------------------------------------------------
    #   Compare reference and fast results
    # -----------------------------------------------------------------

    @test all(isapprox.(yu_ref,  yu_fast ; rtol=1e-6, atol=1e-10))
    @test all(isapprox.(ah_ref,  ah_fast ; rtol=1e-6, atol=1e-10))
    @test all(isapprox.(yu_refAA, yuAA_fast ; rtol=1e-6, atol=1e-8))
    @test all(isapprox.(ah_refAA, ahAA_fast ; rtol=1e-6, atol=1e-8))
end


isinbounds(x,y,z, bounds)=(bounds[1][1] ≤ x ≤ bounds[1][2]) && (bounds[2][1] ≤ y ≤ bounds[2][2]) &&(bounds[3][1] ≤ z ≤ bounds[3][2])

# =============================================================================
#   2️⃣  MULTI‑STEP + BOUNDS INTRA‑CHAIN TEST
# =============================================================================
@testset "IntraChainEnergy – multiple steps + bounds" begin
    # -----------------------------------------------------------------
    #   A simulation with three chains, three time‑steps.
    #   Only *one* of the chains (chain 1) belongs to body == -1,
    #   the other two belong to body = 0 and should therefore be ignored.
    # -----------------------------------------------------------------
    L = 100.0
    Sim = HPSAnalysis.SimData()
    Sim.NAtoms      = 12
    Sim.NChains     = 3
    Sim.NSteps      = 3                                 # three frames
    Sim.ChainStart  = [1,5,9]
    Sim.ChainStop   = [4,8,12]
    Sim.ChainLength = [4,4,4]
    Sim.x           = zeros(12, 3)
    Sim.y           = zeros(12, 3)
    Sim.z           = zeros(12, 3)
    Sim.BoxLength   = [L, L, L]
    Sim.Charges     = zeros(12)
    Sim.IDs         = [1,2,3,4, 1,2,3,4, 1,2,3,4]            # three types, repeated per chain
    Sim.IDToResName = Dict(i=>string('A'+i-1) for i in 1:4)
    Sim.IDToSigmas  = [5,6,7,8]
    Sim.IDToLambdas = collect(LinRange(0.2,1.0,4))
    Sim.Masses      = ones(12)
    Sim.ChainMasses = [4,4,4]
    Sim.EquilibrationTime = 1
    Sim.RGMeasureStep     = 1

    # -----------------------------------------------------------------
    #   Charge pattern (same for all steps)
    # -----------------------------------------------------------------
    id_to_charge = [ +1.0, 0.0, -1.0, 0.5 ]           # A = +1, B = 0, C = –1
    Sim.Charges = [ id_to_charge[id] for id in Sim.IDs ]

    # -----------------------------------------------------------------
    #   Put the three steps at slightly different locations
    # -----------------------------------------------------------------
    for step in 1:3
        # chain 1 – will stay *inside* the bound in every step
        Sim.x[1:4, step] = [30, 34, 38, 42
        ] .+ 0.0
        Sim.y[1:4, step] .= 50
        Sim.z[1:4, step] .= 50

        # chain 2 – moves partially outside the bound in step 2
        Sim.x[5:8, step] = [30, 40, 50, 60] .+ (step == 2 ? -35 : 0)
        Sim.y[5:8, step] .= 50
        Sim.z[5:8, step] .= 50

        # chain 3 – moves partially outside the bound in step 3
        Sim.x[9:12, step] = [30, 40, 50, 60] .+ (step == 3 ? +35 : 0)
        Sim.y[9:12, step] .= 50
        Sim.z[9:12, step] .= 50
    end

    # centre box (so that the minimum‑image convention works)
    Sim.x .-= L/2; Sim.y .-= L/2; Sim.z .-= L/2


    # copy to unwrapped coordinates
    Sim.x_uw= Sim.x
    Sim.y_uw= Sim.y
    Sim.z_uw= Sim.z

    # -----------------------------------------------------------------
    #   Body vector: only chain 1 belongs to body == -1
    # -----------------------------------------------------------------
    #   indices 1‑3 → -1 (active)
    #   indices 4‑9 → 0  (inactive)
    _test_body_vec[] = vcat(fill(-1,4), fill(0,8))

    # -----------------------------------------------------------------
    #   Bounding‑box that *excludes* chain 2 in step 2 and chain 3 in step 3.
    #   The format matches the one used by `computeMeanPairEnergyMatrix_naiv`.
    # -----------------------------------------------------------------
    bounds = [
        (-40.0f0, 40.0f0),   # x ∈ [‑40, 40] Å
        (-40.0f0, 40.0f0),   # y ∈ [‑40, 40] Å
        (-40.0f0, 40.0f0)    # z ∈ [‑40, 40] Å
    ]

    # -----------------------------------------------------------------
    #   Manual reference calculation over *valid* steps only.
    # -----------------------------------------------------------------
    chain_len = Sim.ChainLength[1]             # all chains have same length = 3
    NTypes    = length(keys(Sim.IDToResName))

    yu_ref   = zeros(Float64, chain_len, chain_len)
    ah_ref   = zeros(Float64, chain_len, chain_len)
    yu_refAA = zeros(Float64, NTypes, NTypes)
    ah_refAA = zeros(Float64, NTypes, NTypes)
    n_valid_chains = 0

    Sim.ClusterRange =1:Sim.EquilibrationTime:Sim.NSteps
    for step in Sim.EquilibrationTime : Sim.NSteps
        # -----------------------------------------------------------------
        #   Determine which chains are *completely* inside the user bounds.
        #   Only chains that fulfil the condition for **all** of their atoms
        #   contribute to this time‑step.
        # -----------------------------------------------------------------
        valid_chains = Int[]
        for chain in 1:Sim.NChains
            start_i = Sim.ChainStart[chain]
            stop_i  = Sim.ChainStop[chain]
           
            all_inside= all(map(x->isinbounds(x[1],x[2],x[3], bounds) ,
                            zip(Sim.x[start_i:stop_i, step],
                                Sim.y[start_i:stop_i, step],
                                Sim.z[start_i:stop_i, step]) ))

            if all_inside
                push!(valid_chains, chain)
            end
        end

        # If no chain is valid for this step we simply skip it.
        length(valid_chains) == 0 && continue

        # -------------------------------------------------------------
        #   Accumulate intra‑chain energies **only** for chains that are
        #   inside the bounds *and* belong to body == -1.
        # -------------------------------------------------------------
        for chain in valid_chains
            n_valid_chains += 1
            for (i_rel,i) in enumerate(Sim.ChainStart[chain]:Sim.ChainStop[chain])
                if _test_body_vec[][i]!=-1 continue end
                for (j_rel,j) in enumerate(Sim.ChainStart[chain]:Sim.ChainStop[chain])
                    if !(j_rel > i_rel+1)   continue end
                    if _test_body_vec[][j]!=-1 continue end

                    # distance with minimum‑image convention
                    dx = min_image(Sim.x[i,step] - Sim.x[j,step], L)
                    dy = min_image(Sim.y[i,step] - Sim.y[j,step], L)
                    dz = min_image(Sim.z[i,step] - Sim.z[j,step], L)
                    r  = sqrt(dx*dx + dy*dy + dz*dz)

                    σ = (Sim.IDToSigmas[  Sim.IDs[i] ] + Sim.IDToSigmas[  Sim.IDs[j] ]) / 2.0
                    λ = (Sim.IDToLambdas[ Sim.IDs[i] ] + Sim.IDToLambdas[ Sim.IDs[j] ]) / 2.0

                    en_yu = r < rcut_yu ? yu(r, Sim.Charges[i], Sim.Charges[j]) : 0.0
                    en_ah = r < rcut_ah ? ah(r, σ, λ) : 0.0

                    yu_ref[i_rel, j_rel]   += en_yu
                    ah_ref[i_rel, j_rel]   += en_ah
                    yu_refAA[ Sim.IDs[i], Sim.IDs[j] ] += en_yu
                    ah_refAA[ Sim.IDs[i], Sim.IDs[j] ] += en_ah
                end
            end
        end
    end

    # -----------------------------------------------------------------
    #   Normalise by the number of *valid* frames and symmetrise.
    # -----------------------------------------------------------------
    if n_valid_chains > 0
        yu_ref   ./= n_valid_chains
        ah_ref   ./= n_valid_chains
        yu_refAA ./= n_valid_chains
        ah_refAA ./= n_valid_chains
    end
    symmetrise!(yu_ref);   symmetrise!(ah_ref)
    symmetrise!(yu_refAA); symmetrise!(ah_refAA)

    # -----------------------------------------------------------------
    #   Run the library routine on the same data
    # -----------------------------------------------------------------
    # (pre‑computations that the fast routine expects)
    HPSAnalysis.computeCOM!(Sim)
    HPSAnalysis.computeClustersByChainCOM(Sim; Cutoff=200.0)
    yu_fast, yuAA_fast, ah_fast, ahAA_fast = HPSAnalysis.getIntraChainEnergyMap(Sim, bounds; ϵ_ah=Float32(ϵ_ah))

    # -----------------------------------------------------------------
    #   Final assertions
    # -----------------------------------------------------------------
    @test all(isapprox.(yu_ref  ,   yu_fast ; rtol=1e-6, atol=1e-10))
    @test all(isapprox.(ah_ref  ,   ah_fast ; rtol=1e-6, atol=1e-10))
    @test all(isapprox.(yu_refAA, yuAA_fast ; rtol=1e-6, atol=1e-4 ))
    @test all(isapprox.(ah_refAA, ahAA_fast ; rtol=1e-6, atol=1e-4 ))
    @test n_valid_chains > 0   # sanity‑check that the bounds really filtered something
end