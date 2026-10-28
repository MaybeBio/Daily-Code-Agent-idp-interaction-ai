# Droplet: Agentic Pipeline Design Session

## Executive Summary

**Project**: PhD research on sensing capabilities of biomolecular condensates in gene regulation  
**Current State**: Manual experimental workflow with specification files, sweep scripts, and parallel execution  
**Goal**: Design an AI-agent-driven pipeline for autonomous experimental design, execution, and analysis  
**Constraint**: Maximize use of traditional tools; minimize agent intervention to intelligent tasks only

---

## 1. Scientific Context

### Research Objective
Investigate whether condensate formation and polymer dynamics enable promoters to sense protein concentrations more sensitively than diffusion alone, using lattice-gas simulations of polymer-protein interactions.

### Key Metrics
- Transcription process statistics as a function of concentration
- Fisher information rate (physical limits of sensing)
- Binding and transcription dynamics

---

## 2. Current Workflow Architecture

### 2.1 Specification Files (`spec/`)
**Purpose**: Define complete experiment configuration

**Structure**:
1. **Physical Elements**: Instantiate lattice and species (no parameters)
2. **Simulation Parameters**: 
   - Runtime, burn-in, sampling intervals
   - Parameter sweep definitions
   - Replica count
3. **Target Function**: 
   - Parameter assignment for each sweep point
   - Main Monte Carlo loop
   - Data collection logic
   - Returns results dictionary

### 2.2 Sweep Script (`sim/sweep.jl`)
**Purpose**: Execute parameter sweeps in parallel

**Arguments**:
1. Specification file name
2. CPU key
3. Total CPU count
4. Output directory path

**Behavior**: Loads spec via `include()`, runs sweep loop, writes results

### 2.3 Orchestration (`runJL` bash script)
**Purpose**: Coordinate parallel execution via Unix pipes

### 2.4 Data Organization
**Location**: `data/YYYYMMDD/spec_name_parameters/`

**Example**:
```bash
./runJL sim/sweep.jl spec/prewet_transcription.jl 16 \
  data/20251218/prewet_transcription_32_1300-4000
```

### 3.5 Data visualization
**Location**: viz.ipynb

**Purpose**: Plot and explore simulation data 

---

## 3. Proposed Agentic Pipeline

### 3.1 Phase 1: Experimental Design (CURRENT FOCUS)

**Input**: Natural language experiment description from researcher  
**Output**: Validated experiment proposal

**Required Components**:
- Specification file (`.jl`)
- Experiment summary
- Plot descriptions with expected results
- Technical specifications:
  - Core count
  - Parameter ranges
  - Output data types

**Workflow**:
1. Agent generates proposal from description
2. Iterative refinement with researcher
3. Researcher approval gates execution

### 3.2 Phase 2: Execution (FUTURE)

**Input**: Validated proposal + spec file  
**Output**: Organized simulation results + validation report

**Tasks**:
- Run simulations
- Retrieve and organize data
- Validate output integrity
- Send notifications (start, finish, errors)

### 3.3 Phase 3: Analysis & Reporting (FUTURE)

**Input**: Simulation output data  
**Output**: Polished report with visualizations

**Tasks**:
- Generate plots per proposal
- Create detailed report
- On failure: diagnose issues, propose fixes
- Send completion notification

---

## 4. Design Constraints

### 4.1 Minimize Agent Overhead
- Use traditional software tools wherever possible
- Reserve agent intervention for:
  - Intelligent analysis
  - Hypothesis generation/testing
  - Summarization
  - Debugging
- Automate everything else with scripts

### 4.2 Preserve Repository Structure
- Minimal scaffolding/context additions
- Maintain human readability
- Keep existing file structures intact
- All changes must be understandable by primary developer

### 4.3 Notification System
Implement automated notifications for:
- Simulation start/completion
- Data retrieval completion
- Report readiness
- Error conditions

---

## 5. Agent Task: Initial Assessment

### Your Mission
1. **Explore**: Review the repository structure
2. **Clarify**: Ask questions about unclear aspects
3. **Propose**: Design context and scaffolding for Phase 1 (Experimental Design)
   - How should agents parse experiment descriptions?
   - What templates/schemas are needed for spec file generation?
   - What validation checks should be automated vs. agent-driven?
   - What context files/documentation would help agents understand the codebase?

### Key Questions to Address
- What minimal context do agents need to generate valid spec files?
- How should experiment proposals be formatted for efficient iteration?
- What guardrails prevent agents from breaking existing functionality?
- Where should agent-facing documentation live in the repo?

---

## 6. Success Criteria

A successful Phase 1 design enables:
- ✅ Agent generates syntactically correct spec files from descriptions
- ✅ Proposals contain all required technical details
- ✅ Iteration cycle is efficient (< 3 rounds typical)
- ✅ No disruption to manual workflow capability
- ✅ Repository remains clean and human-maintainable
