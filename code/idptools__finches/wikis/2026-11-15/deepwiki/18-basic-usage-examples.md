# Basic Usage Examples

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [demo/docs_demo/ADBD1.pdb](demo/docs_demo/ADBD1.pdb)
- [demo/docs_demo/epsilon_docs.ipynb](demo/docs_demo/epsilon_docs.ipynb)
- [demo/overview_uses/basic_uses.ipynb](demo/overview_uses/basic_uses.ipynb)

</details>



This page provides simple, practical examples for the most common FINCHES tasks: calculating epsilon values between IDR sequences and generating basic interaction maps. These examples demonstrate the essential workflows using the frontend interfaces with minimal setup.

For comprehensive API documentation, see [Frontend Classes](#5.1). For advanced analysis techniques including phase diagrams and protein-nucleic acid interactions, see [Advanced Analysis Examples](#6.2).

## Setting Up Frontend Interfaces

FINCHES provides two main frontend interfaces that encapsulate different physics models for analyzing IDR interactions. The frontend objects handle all the underlying computational complexity and provide simple methods for common analyses.

```mermaid
graph TB
    User["User Input<br/>Sequences + Parameters"] --> Frontend["Frontend Interface"]
    Frontend --> Mpipi["Mpipi_frontend<br/>Wang-Frenkel + Coulombic"]
    Frontend --> CALVADOS["CALVADOS_frontend<br/>Ashbaugh-Hatch + Yukawa"]
    
    Mpipi --> MpipiParams["salt: float<br/>dielectric: float"]
    CALVADOS --> CalvadosParams["salt: float<br/>pH: float<br/>temp: float"]
    
    MpipiParams --> Results["Epsilon Values<br/>Interaction Maps<br/>Phase Diagrams"]
    CalvadosParams --> Results
```

### Basic Initialization

The simplest way to get started is to initialize the frontend objects with default parameters:

**Default initialization:**
```python
from finches.frontend.mpipi_frontend import Mpipi_frontend
from finches.frontend.calvados_frontend import CALVADOS_frontend

# Use default parameters
mf = Mpipi_frontend()
cf = CALVADOS_frontend()
```

**Custom parameters for physiological conditions:**
```python
# Roughly physiological conditions
mf = Mpipi_frontend(salt=0.150, dielectric=80.0)  # 150 mM salt
cf = CALVADOS_frontend(salt=0.150, pH=7.4, temp=288)  # Physiological pH and ~room temp
```

Sources: [demo/overview_uses/basic_uses.ipynb:95-106](), [demo/docs_demo/epsilon_docs.ipynb:26-32]()

## Basic Epsilon Calculations

Epsilon values provide a scalar measure of the overall attractive or repulsive interaction between two IDR sequences. Lower (more negative) values indicate stronger attractive interactions.

```mermaid
graph LR
    Seq1["IDR Sequence 1"] --> EpsilonCalc["epsilon() method"]
    Seq2["IDR Sequence 2"] --> EpsilonCalc
    EpsilonCalc --> EpsilonValue["Scalar Epsilon Value<br/>(negative = attractive)"]
    
    EpsilonCalc --> Baseline["GS Linker Baseline<br/>(~0.0)"]
    EpsilonCalc --> Corrections["Aliphatic Correction<br/>Charge Correction"]
```

### Single Epsilon Calculation

```python
# Define IDR sequences
med1_idr = 'MESNQSNNGGSGNAALNRGGRYVPPHLRGGDGGAAAAASAGGDDRRGGAGGGGYRRGGGNSGGGGGGGYDRGYNDNRDDRDNRGGSGGYGRDRNYEDRGYNGGGGGGGNRGYNNNRGGGGGGYNRQDRGDGGSSNFSRGGYNNRDEGSDNRGSGRSYNNDRRDNGGD'
med14_idr = 'QDARRRSVNEDDNPPSPIGGDMMDSLISQLQPPPQQQPFPKQPGTSGAYPLTSPPTSYHSTVNQSPSMMHTQSPGNLHAASSPSGALRAPSPASFVPTPPPSSHGISIGPGASFASPHGTLDPSSPYTMVSPSGRAGNWPGSPQVSGPSPAARMPGMSPANPSLHSPVPDASHSPRAGTSSQTMPTNMPPPRKLPQRSWAAS'

# Calculate epsilon value
epsilon_value = mf.epsilon(med1_idr, med14_idr)
print(f"Epsilon value: {epsilon_value}")
```

### Comparing Multiple Interactions

```python
# Compare different IDR pairs
med15_idr = 'KSQASVSDPMNALQSLTGGPAAGAAGIGMPPRGPGQSLGGMGSLGAMGQPMSLSGQPPPGTSGMAPHSMAVVSTATPQTQLQLQQVALQQQQQQQQFQQQQQAALQQQQQQQQQQQFQAQQSAMQQQFQAVVQQQQQLQQQQQQQQHLIKLHHQNQQQIQQQQQQLQRIAQLQLQQQQQQQQQQQQQQQQALQAQPPIQQPPMQQPQPPPSQALPQQLQQMHHTQHHQPPPQPQQPPVAQNQPSQLPPQSQTQPLVSQAQALPGQMLYTQPPLKFVRAPMVVQQPPVQPQVQQQQTAVQTAQAAQMVAPGVQMITEALAQGGMHIRARFPPTTAVSAIPSSSIPLGRQPMAQVSQSSLPMLSSPSPGQQVQTPQSMPPPPQPSPQPGQPSSQPNSNVSSGPAPSPSSFLPSPSPQPSQSPVTARTPQNFSVPSPGPLNTPVNPSSVMSPAGSSQA'

# Calculate multiple interactions
interactions = {
    'med1 + med14': mf.epsilon(med1_idr, med14_idr),
    'med1 + med15': mf.epsilon(med1_idr, med15_idr)
}

for pair, epsilon in interactions.items():
    print(f"{pair}: {epsilon:.3f}")
```

### Understanding Epsilon Order Dependence

Epsilon calculations are "extrinsic" - the order matters because you're effectively "bathing" the first sequence in the second sequence:

```python
# Order matters for epsilon calculations
epsilon_12 = mf.epsilon(med1_idr, med14_idr)  # med1 bathed in med14
epsilon_21 = mf.epsilon(med14_idr, med1_idr)  # med14 bathed in med1

print(f"med1 bathed in med14: {epsilon_12:.3f}")
print(f"med14 bathed in med1: {epsilon_21:.3f}")

# Convert to intrinsic (order-independent) values
intrinsic_12 = epsilon_12 / len(med1_idr)
intrinsic_21 = epsilon_21 / len(med14_idr)
print(f"Intrinsic values are equal: {intrinsic_12:.6f} == {intrinsic_21:.6f}")
```

Sources: [demo/overview_uses/basic_uses.ipynb:163-174](), [demo/overview_uses/basic_uses.ipynb:227-236](), [demo/overview_uses/basic_uses.ipynb:262-266]()

## Controlling Correction Terms

FINCHES applies built-in corrections for aliphatic and charge clustering effects. These can be disabled if needed:

```python
# Default calculation with all corrections
epsilon_default = mf.epsilon(med1_idr, med14_idr)

# Disable specific corrections
epsilon_no_aliphatic = mf.epsilon(med1_idr, med14_idr, use_aliphatic_weighting=False)
epsilon_no_charge = mf.epsilon(med1_idr, med14_idr, use_charge_weighting=False)
epsilon_no_corrections = mf.epsilon(med1_idr, med14_idr, 
                                   use_aliphatic_weighting=False, 
                                   use_charge_weighting=False)

print(f"Default: {epsilon_default:.3f}")
print(f"No aliphatic correction: {epsilon_no_aliphatic:.3f}")
print(f"No charge correction: {epsilon_no_charge:.3f}")
print(f"No corrections: {epsilon_no_corrections:.3f}")
```

Sources: [demo/overview_uses/basic_uses.ipynb:295-309]()

## Epsilon Vector Calculations

While epsilon values provide overall interaction strength, epsilon vectors give residue-level information about attractive and repulsive contributions along the sequence.

```mermaid
graph TB
    SeqPair["IDR Sequence Pair"] --> EpsilonVectors["epsilon_vectors() method"]
    EpsilonVectors --> Attractive["Attractive Vector<br/>numpy array"]
    EpsilonVectors --> Repulsive["Repulsive Vector<br/>numpy array"]
    
    Attractive --> PerResidue["Per-residue<br/>attractive contributions"]
    Repulsive --> PerResidue2["Per-residue<br/>repulsive contributions"]
```

### Basic Vector Calculation

```python
# Calculate epsilon vectors
vectors = mf.epsilon_vectors(med1_idr, med14_idr)
attractive_vector = vectors[0]  # First element is attractive contributions
repulsive_vector = vectors[1]   # Second element is repulsive contributions

print(f"Attractive vector shape: {attractive_vector.shape}")
print(f"First 10 attractive values: {attractive_vector[:10]}")
```

### Analyzing Vector Results

```python
import numpy as np

# Analyze the vectors
total_attractive = np.sum(attractive_vector)
total_repulsive = np.sum(repulsive_vector)
net_interaction = total_attractive + total_repulsive

print(f"Total attractive: {total_attractive:.3f}")
print(f"Total repulsive: {total_repulsive:.3f}")
print(f"Net interaction: {net_interaction:.3f}")

# Find most attractive/repulsive positions
most_attractive_pos = np.argmin(attractive_vector)
most_repulsive_pos = np.argmax(repulsive_vector)

print(f"Most attractive position: {most_attractive_pos}")
print(f"Most repulsive position: {most_repulsive_pos}")
```

Sources: [demo/overview_uses/basic_uses.ipynb:595-597]()

## Comparing Forcefields

FINCHES supports multiple physics models. You can easily compare results between different forcefields:

```mermaid
graph TB
    Sequences["IDR Sequences"] --> MpipiFE["Mpipi_frontend"]
    Sequences --> CalvadosFE["CALVADOS_frontend"]
    
    MpipiFE --> MpipiModel["Mpipi Physics Model<br/>Wang-Frenkel + Coulombic"]
    CalvadosFE --> CalvadosModel["CALVADOS Physics Model<br/>Ashbaugh-Hatch + Yukawa"]
    
    MpipiModel --> Results1["Mpipi Epsilon Values"]
    CalvadosModel --> Results2["CALVADOS Epsilon Values"]
    
    Results1 --> Comparison["Compare Results<br/>Across Forcefields"]
    Results2 --> Comparison
```

### Multi-Forcefield Comparison

```python
# Initialize both frontends
mf = Mpipi_frontend()
cf = CALVADOS_frontend()

# Example sequences (LAF-1 domains)
s1 = 'MESNQSNNGGSGNAALNRGGRYVPPHLRGGDGGAAAAASAGGDDRRGGAGGGGYRRGGGNSGGGGGGGYDRGYNDNRDDRDNRGGSGGYGRDRNYEDRGYNGGGGGGGNRGYNNNRGGGGGGYNRQDRGDGGSSNFSRGGYNNRDEGSDNRGSGRSYNNDRRDNGGD'
s2 = 'LEGMSGDMRSGGGYRGRGGRGNGQRFGGRDHRYQGGSGNGGGGNGGGGGFGGGGQRSGGGGGFQSGGGGGRQQQQQQRAQPQQDWWS'

# Compare forcefield predictions
mpipi_epsilon = mf.epsilon(s1, s2)
calvados_epsilon = cf.epsilon(s1, s2)

print(f"Mpipi epsilon: {mpipi_epsilon:.3f}")
print(f"CALVADOS epsilon: {calvados_epsilon:.3f}")
print(f"Difference: {abs(mpipi_epsilon - calvados_epsilon):.3f}")
```

Sources: [demo/docs_demo/epsilon_docs.ipynb:40-44]()

## GS Linker Baseline

FINCHES uses GS linkers as a baseline reference (epsilon ≈ 0). This provides a neutral reference point for interpreting interaction strengths:

```python
# GS linker baseline
GS_LINKER = 20 * 'GS'  # 40 residue GS repeat

# Check baseline interaction
gs_epsilon = mf.epsilon(GS_LINKER, GS_LINKER)
print(f"GS linker self-interaction: {gs_epsilon:.6f}")
print("This should be close to 0")

# Compare IDR to GS baseline
idr_vs_gs = mf.epsilon(med1_idr, GS_LINKER)
print(f"IDR vs GS linker: {idr_vs_gs:.3f}")
```

Sources: [demo/overview_uses/basic_uses.ipynb:122-124](), [demo/overview_uses/basic_uses.ipynb:199-201]()

## Summary Workflow

The basic FINCHES workflow for epsilon calculations follows this pattern:

```mermaid
graph TD
    Start["Start Analysis"] --> Import["Import Frontend Classes"]
    Import --> Initialize["Initialize Frontend Objects<br/>mf = Mpipi_frontend()<br/>cf = CALVADOS_frontend()"]
    Initialize --> DefineSeqs["Define IDR Sequences<br/>seq1, seq2"]
    DefineSeqs --> CalcEpsilon["Calculate Epsilon<br/>eps = mf.epsilon(seq1, seq2)"]
    CalcEpsilon --> Interpret["Interpret Results<br/>Negative = Attractive<br/>Positive = Repulsive"]
    
    CalcEpsilon --> OptionalVectors["Optional: Calculate Vectors<br/>vectors = mf.epsilon_vectors(seq1, seq2)"]
    OptionalVectors --> OptionalAnalysis["Optional: Analyze Per-residue<br/>Contributions"]
```

This provides the foundation for all FINCHES analyses. More complex workflows build upon these basic epsilon calculations to generate interaction maps, phase diagrams, and other advanced visualizations.

Sources: [demo/overview_uses/basic_uses.ipynb:95-174](), [demo/docs_demo/epsilon_docs.ipynb:26-44]()

---