# Agentic Pipeline Documentation Index

**Purpose**: Quick navigation guide for all Phase 1 documentation

---

## 📋 Start Here

**New to this project?** Read in this order:

1. **`prompts/20260123/phase1_assessment.md`** - Detailed design and questions
2. **`context/spec_template_guide.md`** - How to write spec files
3. **`context/domain_knowledge.md`** - Physical and biological context

---

## 📁 Document Organization

### Planning & Design (`prompts/20260123/`)

| File | Purpose | Length | Audience |
|------|---------|--------|----------|
| `instructions_v2.md` | Original task specification | 200 lines | Context |
| `phase1_assessment.md` | Comprehensive design proposal | 400 lines | Researcher review |
| `proposal_template.yaml` | Experiment proposal schema | 150 lines | Reference |

### Agent Context (`context/`)

| File | Purpose | Length | Agent Use |
|------|---------|--------|-----------|
| `spec_template_guide.md` | Spec file syntax and patterns | 500 lines | Code generation |
| `domain_knowledge.md` | Physics and biology concepts | 600 lines | Proposal design |
| `experiment_types.md` | Existing experiment catalog | 400 lines | Pattern matching |
| `validation_checklist.md` | Quality assurance procedures | 500 lines | Automated validation |

---

## 🎯 Use Cases

### "I want to understand the agent pipeline design"
→ Read: `phase1_assessment.md`

### "I need to write a specification file manually"
→ Read: `context/spec_template_guide.md`, refer to `experiment_types.md`

### "I want to understand the physics"
→ Read: `context/domain_knowledge.md`

### "I'm designing a new experiment"
→ Use: `proposal_template.yaml`, refer to `experiment_types.md` for patterns

### "I need to validate an experiment"
→ Follow: `validation_checklist.md`

---

## 🔑 Key Concepts by Document

### `spec_template_guide.md`
- Specification file structure (3 sections)
- Species types: FastBulk, Bulk, Polymer
- Interaction potentials: `intPot!()`
- Monte Carlo loop patterns
- Results export format
- Common experimental patterns

### `domain_knowledge.md`
- Phase separation & condensation
- Prewetting phenomena
- Cooperative binding
- Hamiltonian, acceptance rate, binding metrics
- Dimensional analysis (lattice size, densities, energies)
- Biological interpretation (TF, RNAP, enhancers)

### `experiment_types.md`
- 7 existing experiment types
- Decision tree for pattern selection
- Modification examples
- Quick reference table

### `validation_checklist.md`
- Syntax validation (Julia compilation)
- Physical validity (density, interactions, timescales)
- Quick functionality test (30-second run)
- Scientific validity (researcher review)
- Validation report template

### `phase1_assessment.md`
- Current workflow analysis
- 20+ clarifying questions
- Proposed agent workflow (5 steps)
- Minimal scaffolding requirements
- Success criteria
- Risk mitigation

---

## 🚀 Quick Reference

### Typical Spec File Structure
```julia
# 1. Physical Parameters
l = Lattice(boxlength=32, nspecies=8)
tf = FastBulk(length=4096)
# ... other species ...

# 2. Simulation Parameters
t_max = 10
step = 1000
alphas = [sweep values]
betas = [0.0]
nparams = length(alphas)

# 3. Target Function
function target!(res, args...)
    slist, alpha, beta = args
    # ... assign parameters ...
    intPot!(species, value)
    init!(l, slist)
    # ... Monte Carlo loop ...
    # ... collect results ...
    return res
end
```

### Typical Interaction Strengths
- Nonspecific: -0.1 to -0.5
- Specific binding: -1.0 to -2.0
- Strong binding: -2.0 to -4.0

### Typical Filling Fractions
- Dilute: ϕ < 0.1
- Moderate: ϕ = 0.1-0.4
- Dense: ϕ = 0.4-0.5
- Overcrowded: ϕ > 0.5 (avoid)

### Typical Runtimes
- Quick test: 0.5 min (30 sec)
- Standard: 10 min
- Long: 20+ min (polymers, dense)

---

## 🔄 Workflow Overview

```
Natural Language Description
         ↓
    [Agent Parsing]
         ↓
   YAML Proposal
         ↓
  [Agent Generation]
         ↓
 Julia Spec File
         ↓
[Automated Validation]
  • Syntax check
  • Physical validity
  • Quick test
         ↓
  Validation Report
         ↓
[Researcher Review]
  • Scientific validity
  • Approval/Iteration
         ↓
  Approved Experiment
         ↓
   [Phase 2: Execution]
         ↓
   [Phase 3: Analysis]
```

---

## 📝 Document Status

| Document | Status | Last Updated |
|----------|--------|--------------|
| All context files | ✅ Complete | 2026-01-23 |
| Proposal template | ✅ Complete | 2026-01-23 |
| Assessment | ✅ Complete | 2026-01-23 |
| Validation tools | ⏳ Pending approval | - |
| Proof of concept | ⏳ Next session | - |

---

## 🎓 Learning Path

**For Agents**:
1. Read `spec_template_guide.md` (syntax and structure)
2. Read `domain_knowledge.md` (physics and biology)
3. Study `experiment_types.md` (existing patterns)
4. Use `validation_checklist.md` (quality assurance)
5. Generate proposals using `proposal_template.yaml`

**For Researchers**:
1. Review `phase1_assessment.md` (design details)
2. Approve/modify context documentation
3. Test workflow with proof of concept
4. Provide feedback for refinement

**For Future Developers**:
1. All documentation is in `context/` and `prompts/`
2. No code dependencies outside `.agent/` (future)
3. Plain text formats (Markdown, YAML, Julia)
4. Version controlled with main repository

---

## 🔗 Related Files in Repository

**Existing Core Files** (not modified):
- `src/*.jl` - Core simulation library
- `spec/*.jl` - Example specification files
- `sim/sweep.jl` - Parallel execution script
- `runJL` - Bash orchestration
- `viz.ipynb` - Visualization notebook
- `Project.toml` - Julia dependencies

**Data Organization**:
- `data/YYYYMMDD/experiment_name/` - Simulation outputs
- `data/YYYYMMDD/experiment_name/output_cpu{N}.json` - Results

---

## 💡 Tips

**Writing Proposals**:
- Start with `proposal_template.yaml`
- Match to existing pattern in `experiment_types.md`
- Validate against `validation_checklist.md`

**Understanding Physics**:
- Read relevant sections of `domain_knowledge.md`
- Check typical parameter ranges
- Review biological interpretation

**Debugging Spec Files**:
- Refer to `spec_template_guide.md` Section 6
- Check common syntax errors
- Verify physical validity thresholds

**Designing New Experiments**:
1. Formulate hypothesis
2. Identify similar pattern in `experiment_types.md`
3. Copy and modify that spec file
4. Validate using checklist

---

## 📧 Next Steps

**Researcher Actions**:
1. Review all context documentation
2. Provide feedback on accuracy and completeness
3. Answer clarifying questions in `phase1_assessment.md`
4. Approve proceeding to proof of concept

**Agent Actions** (when approved):
1. Build validation tools (`.agent/validation_tools.jl`)
2. Reverse-engineer existing spec → YAML
3. Test generation: YAML → spec
4. Demonstrate with new experiment

---

## ✅ Completion Checklist

Phase 1 Design:
- ✅ Repository exploration complete
- ✅ Assessment document created
- ✅ Context documentation written (4 files)
- ✅ Proposal template designed
- ✅ Validation workflow defined
- ✅ No existing files modified
- ⏳ Awaiting researcher review

---

**Last Updated**: January 23, 2026  
**Version**: 1.0  
**Status**: Ready for Review
