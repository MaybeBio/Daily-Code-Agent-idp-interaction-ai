## Bulk

● Great, plan locked in: Metropolis-rate KMC, BulkMonomer-only v1, time-calibrated, sim-time sampling + wall-clock cap. Here's the phased plan.
                                                                                                                                                               
  ---                                                                                                         
  Design decisions (locked before coding)                                                                                                                                                              
  - Rate: k_ij = k0 · min(1, exp(-β·ΔH)). Per-species k0 (sets physical time unit; user sets it from D or attempt-frequency).                                  
  - Move encoding for BulkMonomer: each move is (speciesIdx, particleId, dirIdx ∈ 1..6). Slot ID = stable function of these. 6 slots per particle,             
  pre-allocated.                                                                                                                                               
  - Validity: a move with destination occupied (any species) gets rate 0 — kept in the slot, not removed. Move with rate 0 is just absent from the cumulative  
  sum.                                                                                                                                                         
  - Data structures:
    - rates::Vector{Float64} — slot-indexed, length 6·Σ(species.length).                                                                                       
    - fenwick::Vector{Float64} — Fenwick tree over rates (own minimal helpers in gillespie.jl, no new dep).                                                    
    - The lattice itself is the inverted index: given a site x, the species' latS[x] gives the particle at x in O(1). No new index needed.                     
  - Affected-moves rule (single-species BulkMonomer hop p_old → p_new):                                                                                        
    - dirty_sites = NN1(p_old) ∪ NN1(p_new) (≤14 sites).                                                                                                       
    - For each occupied x ∈ dirty_sites: refresh that particle's 6 moves.                                                                                      
    - For each empty x ∈ dirty_sites: for each y ∈ NN(x), if y is occupied and y ∉ dirty_sites, refresh the one move (particle@y, direction y→x).              
    - Multi-species: same, but iterate over all species at each site.                                                                                          
  - Numerical hygiene: full Q rebuild every N_rebuild events (configurable; default e.g. 10⁶) to bound drift.                                                  
  - Sampling: at fixed Δt_sample of simulation time. Wall-clock cap from target! unchanged.                                                                    
                                                                                                                                                               
  ---                                                                                                                                                          
  File layout                                                                                                                                                  
                  
  - src/gillespie.jl (already a stub) — owns KMCSystem, build/update/sample/apply, Fenwick helpers, KMC loop primitives.
  - src/types.jl — add KMCSystem struct and KMCMove if needed (kept minimal; prefer SoA over a struct per move).                                               
  - test/test_kmc_*.jl — one file per phase.                                                                                                                   
  - spec/kmc_diffusion.jl — minimal v1 spec (a single BulkMonomer doing free diffusion; MSD t¹ check).                                                         
  - sim/sweep.jl — untouched; KMC is a different target! body inside the spec, not a new sweep.                                                                
                                                                                                                                                               
  ---                                                                                                                                                          
  Phased plan                                                                                                                                                  
                  
  Phase 1 — scaffolding & rate primitives (no sampling, no time)
                                                                                                                                                               
  Build:                                                                                                                                                       
  - KMCSystem struct: holds slist, rates, fenwick, t_sim::Float64, Q::Float64, beta::Float64, plus per-species k0 (Vector keyed by species index).             
  - slotId(sysIdx, particleId, dir) — pure, inlinable.                                                                                                         
  - computeRate!(sys, slot) — for the (species, particle, dir) decoded from slot: check destination free across all species in slist; if blocked, rate=0; else
  compute ΔH using existing deltaHam!-style logic (we'll likely need a thin local deltaH_kmc helper since deltaHam! is tightly coupled to Metropolis state — to
   be decided: reuse vs minimal duplication; aim for reuse).                                                                                                   
  - buildMoves!(sys) — full sweep, fills rates, builds Fenwick.                                                                                                
                                                                                                                                                               
  Tests (test_kmc_phase1.jl):                                                                                                                                  
  - A 4³ lattice with 1 BulkMonomer, 1 particle, no neighbors → all 6 rates = k0·1.                                                                            
  - 1 particle with 1 occupied NN, intPot=ε, β=1 → moves toward occupant blocked (rate 0), moves away have rate k0·min(1, exp(+βε)) (gain neighbor lost: ΔH=+ε 
  if intPot is binding=negative; sanity-check sign with closed form).                                                                                          
  - Total Q matches sum of rates.                                                                                                                              
                  
  Gate: pass before Phase 2.                                                                                                                                   
                  
  ---                                                                                                                                                          
  Phase 2 — Fenwick sampling
                                                                                                                                                               
  Build:
  - Minimal Fenwick-tree helpers in gillespie.jl: fenwick_update!, fenwick_total, fenwick_find(u) (returns slot whose cumulative sum first exceeds u).         
  - sampleMove(sys) — draws u ~ U(0, Q), returns slot.                                                                                                         
                                                      
  Tests (test_kmc_phase2.jl):                                                                                                                                  
  - Construct rates [1,2,3,4], sample 10⁵ times with fixed seed, χ² test against expected 0.1/0.2/0.3/0.4.                                                     
  - Edge cases: single nonzero rate → always picks it; rate=0 slots never picked.                                                                              
  - fenwick_total ≈ sum(rates) after random updates.                                                                                                           
                                                                                                                                                               
  Gate: pass before Phase 3.                                                                                                                                   
                                                                                                                                                               
  ---                                                                                                                                                          
  Phase 3 — incremental updates
                                                                                                                                                               
  Build:
  - applyMove!(sys, slot) — translate particle, update lattice, update t_sim += -ln(rand())/Q before the rate refresh (Gillespie order).                       
  - refreshAffected!(sys, p_old, p_new) — implements the dirty-site rule above. Each slot whose rate changes: subtract old, recompute, add new,                
  fenwick_update!.
  - Periodic rebuildQ!(sys) for drift correction.                                                                                                              
                  
  Tests (test_kmc_phase3.jl):                                                                                                                                  
  - After every event for 1000 events, compare incremental Q to a from-scratch buildMoves! recomputation. They must agree to ≤1e-12 relative.
  - A "no-op" move (proposed but ΔH refresh is trivial because no neighbors): no rates besides the mover's own should change.                                  
  - Two-particle test: move particle A; verify particle B's 6 rates are correct vs a hand-computed reference for select configs.
                                                                                                                                                               
  Gate: pass before Phase 4. This is the highest-risk phase — most KMC bugs live here.                                                                         
                                                                                                                                                               
  ---                                                                                                                                                          
  Phase 4 — main loop + sim-time sampling                                                                                                                      
                                                                                                                                                               
  Build:
  - kmc_step!(sys) — sample → advance time → applyMove! → refresh.                                                                                             
  - kmc_run!(sys; t_max_wall, t_burn_in, dt_sample, sample_callback):                                                                                          
    - Outer loop bounded by wall-clock (time_ns()).                  
    - Internal t_next_sample accumulator; whenever t_sim ≥ t_next_sample, invoke callback and bump t_next_sample += dt_sample.                                 
    - Burn-in measured in t_sim, not events.                                                                                                                   
  - Output dict additions: t_sim_final, n_events, dt_sample, k0 per species. Existing fields (runtime_min, boxlength, etc.) preserved.                         
                                                                                                                                                               
  Tests (test_kmc_phase4.jl):                                                                                                                                  
  - Equilibrium parity with Metropolis (the gating test). Single BulkMonomer species, 8³ box, modest density and intPot. Run Metropolis to equilibrium; collect
   histogram of ham. Run KMC with same parameters from the same state; collect same histogram. KS-test or moments agree within MC noise (concrete tolerance:   
  |μ_KMC − μ_M| < 3·max(σ_KMC, σ_M)/√Nbatches).                                                                                                                
  - Free-particle MSD: 1 particle, no interactions, large box. MSD vs t_sim should be linear, slope = 2·D where D = k0 · (1²)/6 for 3D simple-cubic hopping    
  with rate k0/6 per direction (i.e. D = k0/6). Tolerance: 5%.                                                                                             
  - Sim-time advancement sanity: Δt = -ln(u)/Q is correct distribution → mean t_sim per event ≈ 1/⟨Q⟩.                                                         
   
  Gate: parity test must pass.                                                                                                                                 
                  
  ---                                                                                                                                                          
  Phase 5 — multi-species cross-interactions
                                                                                                                                                               
  Once Phase 4 holds, multi-species mostly "just works" because deltaH_kmc already iterates intSpecies. The only new bit is that dirty_sites updates may flip
  rates across species (a TF moving changes rates of an RNAP near it).                                                                                         
                  
  Build:                                                                                                                                                       
  - Generalize buildMoves! and refreshAffected! to iterate over all species in slist for both occupancy checks and rate refresh.
  - Confirm intMat cross-terms enter deltaH_kmc correctly.                                                                                                     
   
  Tests (test_kmc_phase5.jl):                                                                                                                                  
  - Two species, attractive cross-interaction, equilibrium contact-frequency matches Metropolis (same parity protocol as Phase 4).
  - Particle of species A moving toward a B: A's rate to that destination is blocked; A's rate parallel to B is enhanced (binding).                            
                  
  Gate: parity test must pass before writing the spec.                                                                                                         
                  
  ---                                                                                                                                                          
  Phase 6 — spec + driver integration
                                                                                                                                                               
  Build:          
  - spec/kmc_diffusion.jl: copied from the closest existing BulkMonomer spec (likely bulk_diffusion_search.jl), target! body swapped to call kmc_run!. Mirror
  the existing output dict shape (species_name, nframes, boxlength, runtime_min, plus the new KMC fields).                                                     
  - Initial spec runs t_max = 1 minute, ncpus = 1 per CRITICAL_PRINCIPLES.                                
                                                                                                                                                               
  Tests:                                                                                                                                                       
  - The spec executes end-to-end via runJL sim/sweep.jl spec/kmc_diffusion.jl 1 data/.../kmc_diffusion_test.
  - Output JSON loads, has expected keys, MSD-from-coords matches the in-memory MSD from Phase 4.                                                              
                                                                                                 
  Gate: actual runJL run completes and produces non-empty data.                                                                                                
                                                                                                                                                               
  ---                                                                                                                                                          
  Phase 7 — global consistency                                                                                                                                 
                              
  - Run the full test suite once (per style_guide.md, only at the end).
  - Add test_kmc_*.jl includes to runtests.jl.                                                                                                                 
  - Confirm no regression in existing Metropolis tests.                                                                                                        
  - If Polymer KMC (v2) is the next step, this is the natural commit/branch boundary.                                                                          
                                                                                                                                                               
  - Free-particle MSD: 1 particle, no interactions, large box. MSD vs t_sim should be linear, slope = 2·D where D = k0 · (1²)/6
  for 3D simple-cubic hopping with rate k0/6 per direction (i.e. D = k0/6). Tolerance: 5%.
  - Sim-time advancement sanity: Δt = -ln(u)/Q is correct distribution → mean t_sim per event ≈ 1/⟨Q⟩.

  Gate: parity test must pass.

  ---
  Phase 5 — multi-species cross-interactions

  Once Phase 4 holds, multi-species mostly "just works" because deltaH_kmc already iterates intSpecies. The only new bit is that
  dirty_sites updates may flip rates across species (a TF moving changes rates of an RNAP near it).

  Build:
  - Generalize buildMoves! and refreshAffected! to iterate over all species in slist for both occupancy checks and rate refresh.
  - Confirm intMat cross-terms enter deltaH_kmc correctly.

  Tests (test_kmc_phase5.jl):
  - Two species, attractive cross-interaction, equilibrium contact-frequency matches Metropolis (same parity protocol as Phase 4).
  - Particle of species A moving toward a B: A's rate to that destination is blocked; A's rate parallel to B is enhanced
  (binding).

  Gate: parity test must pass before writing the spec.

  ---
  Phase 6 — spec + driver integration

  Build:
  - spec/kmc_diffusion.jl: copied from the closest existing BulkMonomer spec (likely bulk_diffusion_search.jl), target! body
  Build:
  - spec/kmc_diffusion.jl: copied from the closest existing BulkMonomer spec (likely bulk_diffusion_search.jl), target! body swapped to call kmc_run!. Mirror
  the existing output dict shape (species_name, nframes, boxlength, runtime_min, plus the new KMC fields).
  - Initial spec runs t_max = 1 minute, ncpus = 1 per CRITICAL_PRINCIPLES.

  Tests:
  - The spec executes end-to-end via runJL sim/sweep.jl spec/kmc_diffusion.jl 1 data/.../kmc_diffusion_test.
  - Output JSON loads, has expected keys, MSD-from-coords matches the in-memory MSD from Phase 4.

  Gate: actual runJL run completes and produces non-empty data.

  ---
  Phase 7 — global consistency

  - Run the full test suite once (per style_guide.md, only at the end).
  - Add test_kmc_*.jl includes to runtests.jl.
  - Confirm no regression in existing Metropolis tests.
  - If Polymer KMC (v2) is the next step, this is the natural commit/branch boundary.


  ## Polymer (uniform)

  ● Here's the detailed plan for v2 (uniform polymer + polymer–bulk interactions). Same structure as the v1 plan: locked decisions up top, then phased
  deliverables with rigorous tests, then open questions.                                  

  ---                                                                                                                                                          
  Locked design decisions
                                                                                                                                                               
  1. Slot encoding (per-species generalization)                                                                                                              
                                                                                                                                                               
  - Add nslots(s::AbstractSpecies) dispatched function.                                                                                                        
    - nslots(b::BulkMonomer) = 6 * b.length                                                                                                                    
    - nslots(p::Polymer) = 4 * (p.maxlength - 1)  ← pre-allocated to maxlength so add/remove never reallocate.                                                 
  - buildOffsets! switches from 6 * s.length to nslots(s).                                                                                                     
  - Polymer slot decode: (speciesIdx, bondIdx ∈ 0..maxlength-2, dir ∈ 1..4).                                                                                   
                                                                                                                                                               
  2. Stable bondIdx mapping (per-polymer state)                                                                                                                
                                                                                                                                                               
  - One KMCPolymerState struct per polymer in slist, kept in sys.polymer_states::IdDict{AbstractPolymer, KMCPolymerState}:                                     
    - bondIdx_of_key::Dict{Int32, Int}  — bead key → 0-based bondIdx                                                                                         
    - key_at_bondIdx::Vector{Int32} of length maxlength — bondIdx → bead key (0 = free)                                                                        
    - free_bondIdxs::Vector{Int} — stack of free bondIdx values                                                                                                
  - Init: assign bondIdx 0..L-1 to existing keys; rest pushed to free list.                                                                                    
  - addBond! extension: pop 2 free bondIdxs, assign to the new keys (p.keys_vec[end-1:end] after the call).                                                    
  - removeBond! extension: free the bondIdxs of p.ki and p.kj.                                                                                                 
                                                                                                                                                               
  3. Polymer slot active-status                                                                                                                                
                                                                                                                                                               
  A polymer slot is active iff:                                                                                                                                
  - bondIdx is currently assigned to a bead key, AND
  - that bead has a successor (p.pos[ki][5] != 0).                                                                                                             
                                                                                                                                                             
  Otherwise rate = 0. The chain tail's 4 slots are always inactive.                                                                                            
                                                                                                                                                               
  4. Move type from geometry (existing convention)                                                                                                             
                                                                                                                                                               
  Mirroring moveTypes! and proposeMove!:                                                                                                                       
  - dir ∈ 1..4 indexes the orthogonal directions to bondVec.                                                                                                 
  - Compute newPosi, newPosj. Then:                                                                                                                            
    - moveType = 1 (add) iff newPosi != bprev && newPosj != bnext                                                                                            
    - moveType = 2 (kink) iff exactly one equality holds (XOR)                                                                                                 
    - moveType = 3 (remove) iff both hold                                                                                                                      
  - Add Validity: length < maxlength - 1 AND both newPos sites empty across all species.                                                                       
  - Kink validity: the moving bead's new site is empty.                                                                                                        
  - Remove validity: length > 2.                                                                                                                               
                                                                                                                                                               
  Inactive / blocked / out-of-bounds → rate 0.                                                                                                                 
                                                                                                                                                               
  5. k0 field on Polymer                                                                                                                                       
                                                                                                                                                               
  Add k0::Float64 = 1.0 to the Polymer struct (matching BulkMonomer convention). Sets the polymer's bond-move attempt rate.                                    
                                                                                                                                                             
  6. Pure deltaH_kmc(p::Polymer, sys, bondIdx, dir)                                                                                                            
                                                                                                                                                             
  The existing deltaHam!(p::Polymer) mutates p.ham, p.deltaH, p.deltaHP, and partners' deltaHP, so it's not reusable inside KMC's rate-recompute loop. Write a 
  pure version specialized for uniform polymers:                                                                                                             
  - single monomer type, p.patternTails === nothing                                                                                                            
  - assert these in the function body so heterogeneous polymers fail loudly                                                                                    
  - no calls to _addSwapEnergy! / _removeSwapEnergy!                       
  - ~60 lines (3 move types × ~20 lines) of straight ΔH formulas mirroring hamiltonian.jl:258-421 for uniform-only branches                                    
  - returns Inf for any blocked / inactive / out-of-bounds case (callers convert to rate 0)                                                                    
                                                                                                                                                               
  7. Heterogeneous polymers and three-way couplings: out of scope                                                                                              
                                                                                                                                                               
  - Assert p.nMonomerTypes == 1 (or p.pattern == :uniform) in computeRate! and applyMove! polymer branches.                                                    
  - Assert p.monomers[1].couplingMat === nothing (no pattern-tail logic, no on-site couplings).                                                                
                                                                                                                                                               
  8. Refresh strategy for polymer slots                                                                                                                        
                                                                                                                                                               
  Two-tier:                                                                                                                                                    
  - v2.3 simple: refresh all polymer slots after every move (correctness, simplicity). For maxlength ≤ 100, ≤400 slot recomputes per move — acceptable.      
  - v2.3 optimized (deferred unless perf forces it): refresh only polymer slots whose bi, bj, newPosi, or newPosj lies within distance 2 of any dirty site. Add
   _addPolymerSiteAffected!.                                                                                                                                   
                                                                                                                                                               
  For v2.2 (correctness pass), call rebuildQ! after every event to skip refresh entirely. Switch to incremental in v2.3.                                     
                                                                                                                                                               
  9. Geometry helper                                                                                                                                         
                                                                                                                                                               
  The existing identifyBond!(p) requires p.moveIdx indexing into p.keys_vec, but the keys_vec position of a given bead changes under removeBond!. So write a   
  small _polyBondGeometry!(p, ki::Int32) that does the same work as identifyBond! but takes the bead key directly. Two-line difference, no other refactor.
                                                                                                                                                               
  ---                                                                                                                                                        
  Phasing
         
  v2.1 — Polymer slot encoding & rate computation (no moves)
                                                                                                                                                               
  Build
  - nslots(s) dispatched + update buildOffsets!.                                                                                                               
  - KMCPolymerState struct + KMCSystem.polymer_states::IdDict.                                                                                                 
  - _initPolymerState!(sys, p) called from buildMoves! for each polymer in slist.
  - _polyBondGeometry!(p, ki) and _polyMoveVec!(p, dir) helpers.                                                                                               
  - deltaH_kmc(p::Polymer, sys, bondIdx, dir) -> Float64 (pure).                                                                                               
  - computeRate! dispatch: BulkMonomer (existing) vs Polymer (new). Polymer branch decodes slot, calls geometry helpers, validity-checks, returns k0 * min(1,  
  exp(-β·ΔH)) or 0.                                                                                                                                            
  - Add k0 to Polymer struct.                                                                                                                                  
                                                                                                                                                               
  Tests (test_kmc_polymer_phase1.jl)                                                                                                                           
  - Slot count & active-status: line polymer L=8, maxlength=20 → nslots(p) = 76, exactly 28 slots have nonzero rate (4 × 7 active bonds), the chain-tail's 4   
  slots and the 11 unallocated bondIdxs' 44 slots are 0.                                                                                                       
  - Pure-vs-mutating ΔH parity (the gating test): for each (bondIdx, dir):                                                                                     
    a. KMC path: dh_kmc = deltaH_kmc(p, sys, bondIdx, dir)                                                                                                     
    b. Reference path: set p.moveIdx such that keys_vec[moveIdx] = key_at_bondIdx[bondIdx], call identifyBond!, force p.moveType and newPosi/j for the same    
  dir, call deltaHam!(p), capture p.deltaH, then restore p.ham and p.deltaH.                                                                                   
    c. dh_kmc ≈ p.deltaH to floating-point precision.                                                                                                          
  - Hand-computed configs: U-shape, T-shape — verify ΔH for one add, one kink, one remove against pencil-and-paper.                                            
  - Q == sum(rates) and Fenwick consistent with rates after buildMoves!.                                                                                       
                                                                                                                                                               
  Gate before v2.2.                                                                                                                                            
                                                                                                                                                               
  v2.2 — applyMove! for Polymer + bondIdx maintenance                                                                                                          
                                                                                                                                                             
  Build                                                                                                                                                        
  - Dispatch applyMove!(sys, slot):                                                                                                                          
    - Decode → species type → branch.                                                                                                                          
    - Polymer branch: capture KEYCOUNT[] snapshot if moveType==1, run _polyBondGeometry!, _polyMoveVec!, set p.moveType, advance t_sim, call makeMove!(p),
  update bondIdx mapping (push to free for remove, pop for add), update p.length = length(keys(p.pos)).                                                        
    - For now: rebuildQ!(sys) after each move. No incremental refresh.                                                                                         
                                                                                                                                                               
  Tests (test_kmc_polymer_phase2.jl)                                                                                                                           
  - State integrity over event sequences: 1000 events on a polymer + maybe 1 weakly-interacting bulk species. After each event:                                
    - p.length == length(keys(p.pos))                                                                                                                          
    - chain is well-formed: walk from p.head via [5] pointer for p.length steps, ends at 0; back-pointers [4] consistent.                                      
    - Every active polymer slot rate matches a from-scratch computeRate! recomputation.                                                                        
    - Every bulk slot likewise.                                                                                                                                
    - p.latS consistent with p.pos.                                                                                                                            
  - Round-trip: apply an addBond then check that the inverse removeBond slot exists (on the new bond) and would restore the original geometry.                 
  - bondIdx accounting: after N add/remove events, length(state.free_bondIdxs) + length(state.bondIdx_of_key) == maxlength.                                    
                                                                                                                                                               
  Gate before v2.3.                                                                                                                                            
                                                                                                                                                               
  v2.3 — Incremental refresh generalization                                                                                                                    
                                                                                                                                                               
  Build                                                                                                                                                      
  - Generalize _addSiteAffected!:
    - For each site x in dirty set:                                                                                                                            
        - For each bulk species s: existing source-at-x and dest-at-x logic.
      - For each polymer p: enumerate polymer beads at sites within distance 2 of x; for each such bead with an active bondIdx, push all 4 slots to refresh    
  set.                                                                                                                                                         
  - Generalize the dirty-site set:                                                                                                                             
    - Bulk move: existing NN1(p_old) ∪ NN1(p_new).                                                                                                             
    - Polymer add: NN1(newPosi) ∪ NN1(newPosj).                                                                                                                
    - Polymer kink: NN1(old_pos) ∪ NN1(new_pos) of the moving bead.                                                                                            
    - Polymer remove: NN1(bi) ∪ NN1(bj).                                                                                                                       
  - v2.3-simple fallback path (a flag): instead of fine-grained polymer refresh, refresh ALL polymer slots. Useful as a correctness reference and for small    
  polymers.                                                                                                                                                  
                                                                                                                                                               
  Tests (test_kmc_polymer_phase3.jl)                                                                                                                         
  - Incremental == rebuild over 1000 events in three configs:                                                                                                  
    a. Polymer-only, all moves.                                                                                                                                
    b. Polymer + 1 BulkMonomer with self interactions only.                                                                                                    
    c. Polymer + 1 BulkMonomer with cross intMat[poly, bulk].                                                                                                  
  - Fine-grained vs simple-fallback equivalence: same trajectory should produce identical rates whether refresh is "all polymer slots" or "only affected       
  polymer slots". This catches under-coverage in the affected-slot enumeration.                                                                                
                                                                                                                                                               
  Gate before v2.4.                                                                                                                                            
                                                                                                                                                               
  v2.4 — Metropolis–KMC equilibrium parity                                                                                                                     
   
  Two parity testsets, mirroring Phase 4/5 of v1:                                                                                                              
                                                                                                                                                             
  Polymer-only                                                                                                                                                 
  - 16³ box, polymer with maxlength=24, intPot=-0.4, chemPot tuned so equilibrium length ≈ 12, β=1. No bulk.                                                 
  - Run Metropolis (move!(p) loop, ≥200k accepted) → record p.ham and p.length distributions.                                                                  
  - Run KMC → same distributions.                                                            
  - Tolerances: |μ_M − μ_K| < 0.20 for ham, |μ_M − μ_K| < 0.5 for mean length, std agreement ≤ 0.30.                                                           
                                                                                                                                                             
  Polymer + bulk                                                                                                                                               
  - Same polymer + 30 BulkMonomer particles, with intPot!(p_monomer, b, -0.5) cross interaction.                                                             
  - Same protocol. Compare p.ham, b.ham, mean p.length, mean cross-bond count.                                                                                 
                                                                                                                                                             
  Gate before v2.5. The polymer+bulk parity test is the v2 equivalent of v1's gating Phase 4/5 test.                                                           
                                                                                                                                                               
  v2.5 — Spec + driver integration                                                                                                                             
                                                                                                                                                               
  Build                                                                                                                                                        
  - spec/kmc_polymer_msd.jl: copy of kmc_diffusion.jl, swap the species block to:                                                                            
    - one Polymer (uniform, maxlength=64, initShape="line", k0=1)                                                                                              
    - one BulkMonomer (length=200ish, k0=1) — gives polymer–bulk interactions for realism
    - cross intPot                                                                                                                                             
  - Callback samples per bead: position, unwrapped position (single-wrap tracking per bead, similar to v1 spec), p.ham, b.ham, p.length.                       
  - Output JSON additions: polymer_positions, polymer_unwrapped, polymer_length_trajectory.                                                                    
  - Run KMC → same distributions.
  - Tolerances: |μ_M − μ_K| < 0.20 for ham, |μ_M − μ_K| < 0.5 for mean length, std agreement ≤ 0.30.

  Polymer + bulk
  - Same polymer + 30 BulkMonomer particles, with intPot!(p_monomer, b, -0.5) cross interaction.
  - Same protocol. Compare p.ham, b.ham, mean p.length, mean cross-bond count.

  Gate before v2.5. The polymer+bulk parity test is the v2 equivalent of v1's gating Phase 4/5 test.

  v2.5 — Spec + driver integration

  Build
  - Run KMC → same distributions.
  - Tolerances: |μ_M − μ_K| < 0.20 for ham, |μ_M − μ_K| < 0.5 for mean length, std agreement ≤ 0.30.

  Polymer + bulk
  - Same polymer + 30 BulkMonomer particles, with intPot!(p_monomer, b, -0.5) cross interaction.
  - Same protocol. Compare p.ham, b.ham, mean p.length, mean cross-bond count.

  Gate before v2.5. The polymer+bulk parity test is the v2 equivalent of v1's gating Phase 4/5 test.

  v2.5 — Spec + driver integration

  Build
  - spec/kmc_polymer_msd.jl: copy of kmc_diffusion.jl, swap the species block to:
    - one Polymer (uniform, maxlength=64, initShape="line", k0=1)
    - one BulkMonomer (length=200ish, k0=1) — gives polymer–bulk interactions for realism
    - cross intPot
  - Callback samples per bead: position, unwrapped position (single-wrap tracking per bead, similar to v1 spec), p.ham, b.ham, p.length.
  - Output JSON additions: polymer_positions, polymer_unwrapped, polymer_length_trajectory.

  Tests
  - Spec runs end-to-end via runJL sim/sweep.jl spec/kmc_polymer_msd.jl 1 data/.../... on the small (t_max=0.02 min, single-CPU) initial config.
  - Output JSON loads with all required keys plus the polymer extras.
  - Rouse-scaling sanity: for a single internal bead, ⟨Δr²(t)⟩ from the unwrapped trajectory should scale ~`t^0.5at intermediate times before crossing over tot^1.0` (free diffusion of COM) at long times. This is the
  user's stated motivation, so the spec exists explicitly to enable this measurement. The Phase 6 gate checks that the JSON contains the data needed for the dashboard to plot this; the physical validation (correct
  Rouse exponent within fitting noise) is dashboard-side and a follow-up.

  Gate before v2.6.

  v2.6 — Full suite + integration

  - Add test_kmc_polymer_phase{1..4}.jl to runtests.jl.
  - Run full suite. Expect 0 regressions.
  - Update KMC v1 assertion messages from "v1" → "v2 (uniform polymer)".

  ---
  Open questions

  1. Min length for remove. Existing Metropolis allows shrinking to length 2 (1 bond). Should KMC enforce length > 2 or length >= 2 for remove validity? Need to confirm what the existing move!(p) actually allows when
   length=2 → length=0 — is that a valid move? I'll grep the code to confirm before coding, but I lean: enforce length > 2 as KMC validity (i.e., the result has at least 2 beads / 1 bond) since a fully-removed
  polymer leaves the system in a degenerate state.
  2. identifyBond! reuse vs fresh helper. Writing _polyBondGeometry!(p, ki) is ~25 duplicated lines. Alternative: refactor existing identifyBond! to take an optional ki argument with default ki =
  p.keys_vec[p.moveIdx], single function with a 1-line guard. Slight risk of touching existing Metropolis behavior. My preference: refactor. OK with you?
  3. Pure polymer ΔH duplication. Mirror deltaHam!(p::Polymer) as deltaH_kmc(p, sys, bondIdx, dir) — uniform-only branches, ~60 lines. Same risk/benefit tradeoff as the BulkMonomer pure-vs-mutating split in v1. We
  agreed on "pure" for v1 — defaulting to pure for v2 unless you'd rather see the alternative.
  4. ΔH cross-coupling between polymer beads of the same polymer. The existing deltaHam!(p::Polymer) has a special branch for sp isa PolyMonomer && sp.poly === p (different intIdx within the same polymer =
  heterogeneous case). For uniform polymers this branch never fires. Confirming we drop it — single-monomer polymers don't have intra-polymer cross terms beyond the self-self bond counting, which is the standard
  self-interaction term.
  5. Polymer move sampling vs Metropolis Hastings factor. Existing Metropolis has PKINK = 0.8 and (N+2)/(2N+2) / N/(2N+2) Hastings factors in proposeMove! to make DB work with non-uniform proposal density. KMC
  dispenses with all of that — every (bondIdx, dir) has a rate determined solely by ΔH, and DB falls out naturally because forward and reverse rates differ by exactly exp(-βΔH). Just confirming you're OK losing the
  PKINK knob in KMC. It was there to speed up Metropolis equilibration; KMC doesn't need it.
  6. Per-polymer vs global k0. v1 has per-species k0. For polymer-bulk, each species has its own — that's the consistent generalization. But within a polymer, all moves (add/kink/remove) get the same k0. Are you OK
  with that, or do you want separate k0_kink etc. — there's a kinetics literature distinction between local moves (kinks) and reptation-style moves (add/remove at ends). For v2 I'd start with one k0 per polymer; we
  can add knobs in v3 if the dynamics need tuning.

  If you confirm 1–6 (or just say "go"), I'll start with v2.1.