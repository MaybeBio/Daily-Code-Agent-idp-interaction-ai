# Critical Principles for Agentic Experimental Design

**Last Updated**: January 24, 2026  
**Purpose**: Essential rules for AI agents generating experiment proposals

---

## 1. Proposal Format

**Format**: Markdown with embedded YAML code blocks

**Structure**:
```markdown
# Experiment Proposal: [name]

**Date**: YYYYMMDD  
**Status**: Draft

---

## 1. Scientific Goal

```yaml
hypothesis: "..."
research_question: "..."
expected_outcome: "..."
```

---

## 2. Physical Setup
... (continue with all 7 sections)
```

**Template**: See `prompts/20260123/proposal_template.md`

**Examples**:
- `prompts/20260123/proposal_prewet_enhancers.md` (reference)
- `prompts/20260123/proposal_poly_bulk_phases.md` (test case)

---

## 2. Start Small, Scale Up

**ALWAYS begin with minimal resources**:

```yaml
runtime:
  t_max: 1  # minute (NOT 5, NOT 10)
  
execution:
  cpu_count: 1  # (NOT 16, NOT parallel)
```

**Why**:
- Fast feedback loop (errors, visualization, physics)
- Iterate on proposal before committing compute time
- Catch bugs early (syntax, initialization, metrics)

**Override rule**: When copying from existing spec file, REPLACE the template's runtime/CPU values with small initial values, regardless of what the original used.

**Scaling sequence**:
1. **Test**: 1 min, 1 CPU → validate proposal/spec/dashboard
2. **Iterate**: Fix errors, refine parameters
3. **Medium**: 5-10 min, 1-4 CPUs → verify equilibration
4. **Production**: 10+ min, 16+ CPUs → final data

---

## 3. Minimize Changes

**When creating new spec files**:
- Start from most similar existing spec file
- Copy-paste liberally
- Modify ONLY what's necessary
- Keep same structure, naming, patterns

**Workflow**:
1. Identify closest template (check `experiment_types.md`)
2. Copy that spec file
3. Make minimal modifications
4. Test immediately

**Don't**: Write spec files from scratch  
**Do**: Copy, modify minimally, test

---

## 4. Workflow Summary

```
Natural Language Description
         ↓
    [Ask clarifying questions]
         ↓
  Generate Proposal (markdown + YAML)
         ↓
  Generate Spec File (copy + modify from template)
         ↓
    [User tests: 1 min, 1 CPU]
         ↓
  Create Dashboard (based on proposal plot specs)
         ↓
  Add entry to viz.ipynb
         ↓
    [Iterate based on results]
         ↓
    [Scale up when validated]
```

---

## 5. File Naming Conventions

**Proposals**:
- Location: `prompts/YYYYMMDD/`
- Format: `proposal_[experiment_name].md`
- Example: `proposal_poly_bulk_phases.md`

**Spec files**:
- Location: `spec/`
- Format: `[experiment_name].jl`
- Example: `poly_bulk_phases.jl`

**Dashboard files**:
- Location: Root directory
- Format: `dashboard_[experiment_name].py`
- Example: `dashboard_poly_bulk_phases.py`

**Data output**:
- Location: `data/YYYYMMDD/`
- Format: `[experiment_name]_[boxlength]_[param_info]/`
- Example: `data/20260124/poly_bulk_phases_16_n1-1/`

---

## 6. Context Files Reference

**Core documentation**:
- `context/spec_template_guide.md` - Spec file syntax reference
- `context/visualization_guide.md` - Dashboard patterns
- `context/domain_knowledge.md` - Physics and biology background
- `context/CRITICAL_PRINCIPLES.md` - This file

**Proposal workflow**:
- `prompts/20260123/proposal_template.md` - Standard template
- `prompts/20260123/proposal_prewet_enhancers.md` - Reference example
- `prompts/20260123/instructions.md` - Original agent instructions

---

## 7. Common Mistakes to Avoid

❌ **Don't**:
- Generate proposals without asking clarifying questions
- Start with long runtimes (10+ min) for initial tests
- Use parallel execution (ncpus > 1) before validating
- Write spec files from scratch
- Create pure YAML proposals (use markdown + YAML blocks)

✅ **Do**:
- Ask questions when description is ambiguous
- Always start: 1 min, 1 CPU
- Copy from similar spec file (minimize changes)
- Use markdown format with embedded YAML
- Test → iterate → scale up

---

**End of Principles**
