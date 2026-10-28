# Phase 1 Assessment: Agentic Experimental Design Pipeline

**Date**: January 23, 2026  
**Branch**: `agentic-exp-design`  
**Status**: Initial Assessment Complete

---

## 1. Repository Exploration Summary

### 1.1 Current Architecture Understanding

**Specification Files** (`spec/*.jl`):
- Self-contained Julia scripts defining complete experiments
- Three main sections:
  1. **Physical Parameters**: Lattice + species instantiation (using constructors like `FastBulk`, `Bulk`, `Polymer`)
  2. **Simulation Parameters**: Runtime, sweep ranges, parallelization setup
  3. **Target Function**: Monte Carlo loop with data collection logic

**Key Patterns Observed**:
```julia
# Physical setup
l = Lattice(boxlength=32, nspecies=8)
tf = FastBulk(length=4096, initShape="uniform")
dna = FastBulk(length=13, initShape="line", initPos=[1,10,16])

# Sweep configuration
vals = range(1300, stop=4000, step=180)
alphas = [repeat each value nrepeats times]
nparams = length(alphas)

# Target function structure
function target!(res, args...)
    # 1. Unpack parameters
    # 2. Set interaction potentials via intPot!()
    # 3. Initialize with init!(l, slist)
    # 4. Monte Carlo loop with move!()
    # 5. Collect metrics (hamiltonian, binding, transcription, density)
    # 6. Return results dict
end
```

**Execution Flow**:
```bash
./runJL sim/sweep.jl spec/prewet_transcription.jl 16 data/20251218/output_dir
# Spawns 16 parallel Julia processes via Unix pipes
# Each CPU handles subset of parameter sweep
# Outputs to data/YYYYMMDD/experiment_name_params/output_cpu{N}.json
```

**Data Organization**:
- Directory structure: `data/YYYYMMDD/spec_name_param_range/`
- Output format: `output_cpu{N}.json` containing nested dicts:
  - Level 1: Parameter sweep index
  - Level 2: Beta value index
  - Level 3: Results dict with metrics

---

## 2. Clarifying Questions

### 2.1 Scientific Context

**Q1**: What is the typical experiment lifecycle?
- How many iterations of refinement are normal before a "good" experiment?
- What are common failure modes (bugs, unexpected results, computational issues)?

**A1**: The experiment lifecycle typically starts during the development of a new feature in `dev.jl`. Once tested, the script used in this early iteration phase typically gets copied and modified in a spec file, since new experiments often require new features to be implemented (this might change as the code becomes more comprehensive). For example, once the multicomponent hamiltonian feature was done, I set up two parallel rods with promoters with stronger interaction potentials in dev.jl, and iterated by running quick simulations and visualizing them in viz.ipynb, in the 20251210 -> Prewet Enhancers section of the notebook. I set up a custom dashboard to visualize the data I'm interested in, in this case the acceptance rate, hamiltonian, and YZ and XZ cross sectional occupation densities of the lattice, to observe the condensation effect as the rods distance was varied. Once this loop was closed, I copy the script from dev.jl to its own specification file called prewet_enhancers.jl, and make small test runs on a single CPU, iterating until errors are gone. When that is done, I will scale up the simulation, adding more range to the parameters, more CPUs, more runtime, and if required, send it to the cluster for overnight crunching.

- From the point when the specification file is created and iteration has moved past the dev stage, I would say anywhere between 5 to 50 iterations are required, sometimes going back to dev to fix or add a feature before getting the spec file right.

- Common failure modes are hamiltonian not converging due to insufficient simulation steps, data read/write errors along the json write (julia side, in the spec file) to json read (python side, in viz.ipynb) sequence, syntax errors in the speficiation file, swept parameter range doesn't showcase the desired output behavior (incorrectly chosen), etc. These are just examples, but failure modes are very diverse and can come from many different places.

**Q2**: What domain knowledge should agents understand?
- Statistical mechanics concepts (Hamiltonians, chemical potentials)?
- Phase transition phenomenology (droplet formation, prewetting)?
- Polymer physics (kink moves, self-avoidance)?

**A2**: Flory-Huggins and Cahn-Hillard theories of phase separation, pre-wetting, terms like 'promoter', 'enhancer', 'transcription factor', 'RNA Polymerase', 'mRNA transcription', and their biological meaning in the context of gene regulation, Metropolis Hastings, general Markov-Chain Monte carlo knowledge, ball and chain models of self-avoiding polymers with kink, add and remove moves.

**Q3**: What are the key experimental design decisions?
- Choosing between species types (`Bulk` vs `FastBulk` vs `Polymer`)?
- Setting interaction potential magnitudes (what's physically meaningful)?
- Determining appropriate sweep ranges and resolutions?
- Selecting runtime/sampling intervals for convergence?

**A3**: 
- Bulk is a legacy type, choose FastBulk by default. Polymer is chosen when a polymer is required in the simulation.
- For now this is arbitrary, but the magnitudes should follow from general biological domain knowleedge and the qualitative description given by the researcher. For example, the interaction between an enhancer and a transcription factor should be large compared to that between a TF and a generic location along the Polymer. TF and RNA Polymerase should interact strongly too. TFs have a self-interaction of a certain strength to model condensation. These are usually parameters set by the researcher explicitly.
- This is set by the researcher explicitly in the experiment description
- Also set by the researcher



### 2.2 Technical Workflow

**Q4**: How are experiments currently initiated?
- Do you write natural language notes first, then code?
- Are there experimental "templates" you reuse?
- How much copy-paste between similar experiments?

**A4**:
- I write notes from scientific meetings, and then translate them to code immediately without going through a structured experimental proposal or design document/template.
- I don't use any template, but I usually use as much as I can from previous speficiation files and modify them as little as possible
- As much as possible

**Q5**: What validation checks catch most errors?
- Syntax errors (Julia compiler)?
- Runtime errors (NaNs, infinite loops)?
- Scientific validity (unphysical parameters, poor sampling)?

**A5**: Yes, syntax errors. Recently I have been setting a maximum wall-clock time as a termination criterion for simulations, so I don't worry about infinite or slow runtimes anymore. The problem is now to make sure equilibration happens in the set time, and this is done by iterating and building empirical knowledge, scaling up simulations progressively.
I usually record a sample (I call it 'frame') after a large number of monte-carlo moves (e.g. 10000), so this is not usually an error I deal with.

**Q6**: What analysis patterns are most common?
- Time series plots (binding/transcription dynamics)?
- Phase diagrams (density vs parameters)?
- Statistical measures (Fisher information, correlations)?

**A6**: I typically make parametric plots of some quantity as a function of a varied parameter. For example in the bulk_phases and poly_phases experiments I plotted density and gyration radius respectively versus temperature. I haven't plotted time-series or statistical measures yet, but this is the objective for upcoming experiments.

### 2.3 Agent Interaction Preferences

**Q7**: How detailed should experiment proposals be?
- Just the spec file + brief summary?
- Full technical writeup with expected results?
- Interactive refinement vs batch approval?

**A7**: Just the spec file, a brief summary of the key elements being tested, and changes with respect to the most similar existing spec file or the previous proposal iteration. I will tell the agent to iterate on the proposal multiple times before approving it, asking for small validation runs. Full approval means the experiment is ready for large-scale deployment on a cluster, and the proposal will have to change many times before that happens.

**Q8**: What level of automation for validation?
- Should agents run syntax checks automatically?
- Quick test runs on small grids to verify code?
- Automated plot generation for proposal visualization?

**A8**: 
- No syntax checks run by agents, let the compiler handle this.
- Yes, quick test runs on small grids to verify code is a good idea
- Maybe, if I ask for it explicitly.

---

## 3. Proposed Phase 1 Design

### 3.1 Context Architecture

**Minimal Agent Context Files** (in `context/`):

1. **`spec_template_guide.md`**
   - Canonical structure of spec files
   - Available species types and their parameters
   - Common interaction potential patterns
   - Typical parameter ranges for different experiment types

2. **`domain_knowledge.md`**
   - Physical system overview (condensates, polymers, lattice gas)
   - Key metrics and their interpretations
   - Dimensional analysis guidelines (choosing boxlength, particle counts)
   - Interaction strength conventions

3. **`experiment_types.md`**
   - Catalog of existing experiments with descriptions
   - Phase diagrams (bulk_phases, poly_phases, flory_huggins)
   - Prewetting studies (prewet_transcription, prewet_enhancers)
   - When to use each pattern

4. **`validation_checklist.md`**
   - Syntax requirements (Julia conventions)
   - Physical validity checks (energy scales, densities)
   - Computational feasibility (runtime estimates)
   - Data format requirements

### 3.2 Experiment Proposal Schema

**Structured Proposal Format** (YAML for human readability):

```yaml
experiment:
  name: "prewet_transcription_sensitivity"
  date: "20260123"
  description: |
    Test how transcription factor concentration affects RNA polymerase 
    binding dynamics at promoter sites with enhancer condensates.

scientific_goal:
  hypothesis: "Higher TF concentration increases RNAP binding events"
  metrics:
    - binding_frequency: "TF-enhancer binding time series"
    - transcription_events: "RNAP-promoter binding time series"
    - density_profile: "TF spatial distribution around enhancers"

physical_setup:
  lattice:
    boxlength: 32
    nspecies: 8
  species:
    - name: "tf"
      type: "FastBulk"
      length: "SWEEP"  # 1300-4000 in 180 steps
      initShape: "uniform"
    - name: "rnap"
      type: "FastBulk"
      length: 32
      initShape: "uniform"
    - name: "dna_elements"
      components: ["dna1", "enhancer", "dna2", "promoter"]
      geometry: "parallel_filaments"
  
  interactions:
    - species: ["tf", "self"]
      potential: -1.0
      rationale: "Attractive for condensation"
    - species: ["tf", "enhancer"]
      potential: -2.0
      rationale: "Strong specific binding"
    - species: ["rnap", "promoter"]
      potential: -2.0
      rationale: "Transcription initiation"
    - species: ["tf", "rnap"]
      potential: -2.0
      rationale: "Cooperative interaction"

simulation_parameters:
  runtime_minutes: 10
  burn_in_fraction: 0.8
  sampling_interval: 1000
  sweep:
    parameter: "tf.length"
    range: [1300, 4000]
    step: 180
    replicates: 1
  
execution:
  cpu_count: 16
  estimated_runtime: "~30 minutes"
  output_directory: "data/20260123/prewet_transcription_32_1300-4000"

expected_results:
  plots:
    - type: "time_series"
      data: ["binding", "transcription"]
      description: "Binding events vs time for each TF concentration"
    - type: "scatter"
      x: "tf_concentration"
      y: "mean_binding_frequency"
      description: "Dose-response curve"
    - type: "heatmap_3d"
      data: "tf_density_lattice"
      description: "TF spatial distribution at highest concentration"
  
  validation:
    - check: "binding_increases_with_tf"
      threshold: "Pearson r > 0.7"
    - check: "no_binding_without_tf"
      threshold: "binding_events < 5 when tf.length < 100"

technical_notes:
  - "Using FastBulk for speed (no polymer constraints)"
  - "Parallel DNA filaments separated by 3 lattice units"
  - "Monitor acceptance rates for equilibration"
```

### 3.3 Agent Workflow (Phase 1)

**Input**: Natural language experiment description from researcher

**Agent Tasks**:

1. **Parsing & Understanding**
   - Extract scientific goal, key species, interaction logic
   - Identify which existing experiment pattern is closest match
   - Map natural language to domain concepts using `domain_knowledge.md`

2. **Proposal Generation**
   - Create YAML proposal using schema
   - Generate spec file (`spec/*.jl`) from proposal
   - Use `spec_template_guide.md` for correct syntax

3. **Automated Validation**
   - **Syntax check**: Run `julia --check-bounds=yes spec/proposal.jl`
   - **Quick test**: Run 1 CPU, 1 parameter, short runtime (30 sec)
   - **Physical sanity**: Verify particle counts < boxlength³, energy scales reasonable
   - Report validation results with proposal

4. **Iterative Refinement**
   - Present proposal + validation results to researcher
   - Incorporate feedback (parameter adjustments, interaction changes, etc.)
   - Regenerate spec file and revalidate
   - Typical cycle: 2-3 iterations

5. **Final Approval Gate**
   - Researcher reviews complete proposal
   - Explicit approval required before Phase 2 (execution)
   - Spec file committed to `spec/` with descriptive name

---

## 4. Minimal Scaffolding Requirements

### 4.1 New Files (All Human-Readable)

```
context/
├── spec_template_guide.md       # How to write spec files
├── domain_knowledge.md           # Physics concepts
├── experiment_types.md           # Existing patterns
└── validation_checklist.md       # Quality checks

prompts/
└── 20260123/
    ├── instructions_v2.md        # (existing)
    ├── phase1_assessment.md      # (this document)
    └── proposal_schema.yaml      # Template for proposals

.agent/
├── proposal_template.yaml        # Blank proposal form
└── validation_tools.jl           # Helper functions for syntax/quick tests
```

### 4.2 No Changes to Existing Files

- Core library (`src/*.jl`): **Untouched**
- Existing specs (`spec/*.jl`): **Preserved**
- Sweep script (`sim/sweep.jl`): **Unchanged**
- Data organization: **Maintained**

### 4.3 Notification Hooks (Future)

Placeholder locations for Phase 2/3:
```bash
.agent/hooks/
├── notify_start.sh      # Email/Slack on experiment start
├── notify_complete.sh   # On successful completion
└── notify_error.sh      # On failures
```

---

## 5. Validation Strategy

### 5.1 Automated Checks (Agent-Driven)

**Syntax Validation**:
```bash
julia --check-bounds=yes --project=@. -e 'include("spec/proposal.jl")'
# Should compile without errors
```

**Quick Functionality Test**:
```julia
# Modify proposal to run minimal case
alphas = [alphas[1]]  # Single parameter
betas = [betas[1]]
t_max = 0.5  # 30 seconds
# Run and verify output JSON structure
```

**Physical Sanity Checks**:
- Particle count ≤ 0.5 × boxlength³ (avoid overcrowding)
- |interaction_potentials| ≤ 5.0 (typical energy scale)
- Runtime > 2 minutes (meaningful statistics)
- Sampling interval < total_steps / 100 (sufficient data points)

### 5.2 Researcher Validation (Required)

**Approval Checklist**:
- [ ] Scientific hypothesis is testable with proposed setup
- [ ] Parameter ranges cover physically interesting regime
- [ ] Interaction potentials reflect intended biology
- [ ] Computational resources reasonable (not weeks of runtime)
- [ ] Output metrics align with analysis goals

**Interactive Review**:
- Agent presents proposal in Markdown + YAML
- Researcher can request parameter tweaks via natural language
- Agent regenerates spec file, reruns validation
- Iterate until satisfactory

---

## 6. Success Criteria (Phase 1)

### 6.1 Functional Requirements

✅ **Agent generates valid spec files**:
- Syntactically correct Julia code
- Follows established patterns from `spec/`
- Compiles without errors

✅ **Proposals capture all necessary details**:
- Complete physical setup (species, interactions)
- Clear simulation parameters (sweep ranges, runtime)
- Expected outputs and validation criteria
- Execution metadata (CPU count, directory)

✅ **Efficient iteration cycle**:
- < 5 minutes per iteration (generation + validation)
- Typical approval in ≤ 3 rounds
- Clear feedback mechanism

✅ **Non-disruptive to manual workflow**:
- Can still write spec files by hand
- No dependencies on agent tools for core functionality
- Agent-generated specs indistinguishable from human-written

✅ **Clean repository**:
- New files limited to `context/`, `prompts/`, `.agent/`
- Documentation is human-maintainable
- No auto-generated cruft

### 6.2 Quality Metrics

- **Syntax error rate**: < 5% on first proposal generation
- **Physical validity**: 100% (automated checks prevent unphysical setups)
- **Researcher satisfaction**: Useful proposals that require only minor tweaks
- **Time savings**: 30-50% reduction in experiment setup time

---

## 7. Next Steps

### 7.1 Immediate Actions (Agent Tasks)

1. **Create context documentation**:
   - [ ] Write `spec_template_guide.md` by analyzing existing specs
   - [ ] Extract domain knowledge into `domain_knowledge.md`
   - [ ] Catalog experiment types in `experiment_types.md`
   - [ ] Draft `validation_checklist.md`

2. **Implement proposal schema**:
   - [ ] Formalize YAML structure in `proposal_schema.yaml`
   - [ ] Create blank template in `.agent/proposal_template.yaml`

3. **Build validation tools**:
   - [ ] Write `.agent/validation_tools.jl` with helper functions
   - [ ] Test syntax checking on existing specs
   - [ ] Prototype quick-test runner

4. **Proof of concept**:
   - [ ] Take an existing spec (e.g., `prewet_transcription.jl`)
   - [ ] Generate proposal YAML from it (reverse engineering)
   - [ ] Regenerate spec file from YAML
   - [ ] Verify equivalence

### 7.2 Researcher Review Points

**Decision Points**:
1. Approve YAML schema structure (too verbose? missing fields?)
2. Review context documentation for accuracy and completeness
3. Test proposal generation with a new experiment idea
4. Iterate on workflow until natural and efficient

**Questions for Next Session**:
- Which existing experiment should be the first test case?
- Are there specific interaction patterns to highlight in templates?
- What level of plot specification is helpful in proposals?

---

## 8. Risks & Mitigations

### 8.1 Identified Risks

**Risk**: Agent generates syntactically valid but scientifically nonsensical experiments
- **Mitigation**: Mandatory researcher approval gate + physical sanity checks

**Risk**: Context documentation becomes stale as codebase evolves
- **Mitigation**: Version context files, include last-updated dates, periodic review

**Risk**: Proposal schema too rigid, doesn't capture experimental creativity
- **Mitigation**: Allow "custom_parameters" section for novel setups, iterative refinement

**Risk**: Validation tools introduce dependencies that complicate manual workflow
- **Mitigation**: Keep validation tools separate in `.agent/`, never import into `src/`

### 8.2 Open Questions

- How to handle truly novel experiments that don't fit existing patterns?
- Should agents suggest "similar experiments" from literature/data?
- How to incorporate failed experiments into agent learning?

---

## 9. Long-Term Vision

### Phase 2 Preview: Execution Automation
- Monitor simulation progress
- Detect hangs/errors and restart
- Validate output data integrity
- Send notifications at key milestones

### Phase 3 Preview: Intelligent Analysis
- Auto-generate plots from proposal specifications
- Statistical analysis (correlations, Fisher information)
- Failure diagnosis (poor sampling, unequilibrated, bugs)
- Report generation with LaTeX-quality figures

### Eventual Integration
- Closed-loop experimentation: agent proposes follow-ups based on results
- Hypothesis refinement: update proposals based on data
- Literature integration: suggest experiments from papers
- Collaboration: multi-agent teams for complex projects

---

**End of Assessment**

**Prepared by**: AI Agent  
**For**: PhD Researcher (Droplet Project)  
**Awaiting**: Feedback and approval to proceed with context documentation
