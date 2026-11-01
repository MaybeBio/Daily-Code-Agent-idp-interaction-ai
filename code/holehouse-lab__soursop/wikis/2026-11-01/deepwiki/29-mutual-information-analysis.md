# Mutual Information Analysis

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [soursop/_internal_data.py](soursop/_internal_data.py)
- [soursop/configs.py](soursop/configs.py)
- [soursop/soursop.py](soursop/soursop.py)
- [soursop/ssexceptions.py](soursop/ssexceptions.py)
- [soursop/ssio.py](soursop/ssio.py)
- [soursop/ssmutualinformation.py](soursop/ssmutualinformation.py)
- [soursop/ssnmr.py](soursop/ssnmr.py)
- [soursop/tests/test_ssmutual_information.py](soursop/tests/test_ssmutual_information.py)

</details>



This page documents the mutual information (MI) calculation functionality provided by the `ssmutualinformation` module. Mutual information is a statistical measure that quantifies the amount of information shared between two variables, making it useful for identifying correlations between different structural observables in molecular dynamics trajectories.

For other specialized analysis capabilities, see [NMR Chemical Shift Prediction](#6.1), [PRE Profile Calculation](#6.2), and [Polymer Physics Utilities](#6.4). For general structural analysis methods, see [SSProtein: Single Protein Analysis](#4).

**Sources:** [soursop/ssmutualinformation.py:14-18]()

---

## Module Architecture

Unlike other specialized analysis modules, `ssmutualinformation` is implemented as a collection of standalone functions rather than a class. This design allows these statistical utilities to be used flexibly with any observables extracted from trajectory analysis.

```mermaid
graph TB
    subgraph SSProtein["SSProtein Analysis Methods"]
        RG["get_radius_of_gyration()"]
        E2E["get_end_to_end_distance()"]
        Distances["get_inter_residue_distances()"]
        SASA["get_SASA()"]
        Other["50+ other methods..."]
    end
    
    subgraph Observables["Observable Arrays"]
        X["X: np.array<br/>Observable 1<br/>(e.g., Rg values)"]
        Y["Y: np.array<br/>Observable 2<br/>(e.g., E2E values)"]
        W["weights: np.array<br/>(optional)"]
    end
    
    subgraph ssmutualinformation["ssmutualinformation.py"]
        CalcMI["calc_MI(X, Y, bins,<br/>weights=False,<br/>normalize=False)"]
        ShanEntropy["shan_entropy(c)"]
    end
    
    subgraph Results["Statistical Measures"]
        MI["Mutual Information:<br/>I(X;Y) = H(X) + H(Y) - H(X,Y)"]
        NMI["Normalized MI:<br/>NMI = I(X;Y) / H(X,Y)"]
    end
    
    RG --> X
    E2E --> Y
    Distances --> X
    SASA --> Y
    Other --> X
    Other --> Y
    
    X --> CalcMI
    Y --> CalcMI
    W --> CalcMI
    
    CalcMI --> ShanEntropy
    ShanEntropy --> MI
    CalcMI --> MI
    CalcMI --> NMI
    
    style ssmutualinformation fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:1-136]()

---

## Information Theory Foundations

Mutual information quantifies the reduction in uncertainty about one variable when the other is known. The calculation is based on Shannon entropy and joint probability distributions.

### Mathematical Framework

| Measure | Formula | Interpretation |
|---------|---------|----------------|
| **Shannon Entropy** | `H(X) = -Σ p(x) log p(x)` | Uncertainty in distribution X |
| **Joint Entropy** | `H(X,Y) = -Σ p(x,y) log p(x,y)` | Uncertainty in joint distribution |
| **Mutual Information** | `I(X;Y) = H(X) + H(Y) - H(X,Y)` | Shared information between X and Y |
| **Normalized MI** | `NMI = I(X;Y) / H(X,Y)` | MI scaled to [0,1] range |

```mermaid
flowchart LR
    subgraph Input["Input Data"]
        X["X observable:<br/>np.array of values"]
        Y["Y observable:<br/>np.array of values"]
        Bins["bins:<br/>histogram binning"]
    end
    
    subgraph Histogramming["Histogram Generation"]
        Hist2D["np.histogram2d(X, Y, bins)<br/>→ c_XY"]
        HistX["np.histogram(X, bins)<br/>→ c_X"]
        HistY["np.histogram(Y, bins)<br/>→ c_Y"]
    end
    
    subgraph Entropy["Entropy Calculation"]
        HX["shan_entropy(c_X)<br/>→ H_X"]
        HY["shan_entropy(c_Y)<br/>→ H_Y"]
        HXY["shan_entropy(c_XY)<br/>→ H_XY"]
    end
    
    subgraph MI["Mutual Information"]
        Calc["MI = H_X + H_Y - H_XY"]
        Norm["NMI = MI / H_XY"]
    end
    
    X --> Hist2D
    Y --> Hist2D
    X --> HistX
    Y --> HistY
    Bins --> Hist2D
    Bins --> HistX
    Bins --> HistY
    
    Hist2D --> HXY
    HistX --> HX
    HistY --> HY
    
    HX --> Calc
    HY --> Calc
    HXY --> Calc
    Calc --> Norm
    HXY --> Norm
    
    style Entropy fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:26-105](), [soursop/ssmutualinformation.py:110-136]()

---

## Function Reference

### calc_MI()

Calculate mutual information between two observables with optional weighting and normalization.

```mermaid
graph TD
    Function["calc_MI(X, Y, bins,<br/>weights=False,<br/>normalize=False)"]
    
    subgraph Validation["Input Validation"]
        LengthCheck["len(X) == len(Y)?"]
        BinMin["min(bins) ≤ min(X, Y)?"]
        BinMax["max(bins) ≥ max(X, Y)?"]
    end
    
    subgraph Processing["Processing Path"]
        Weighted{"weights<br/>provided?"}
        Hist["Generate histograms<br/>with/without weights"]
        Shannon["Calculate Shannon<br/>entropies H_X, H_Y, H_XY"]
        MICalc["MI = H_X + H_Y - H_XY"]
    end
    
    subgraph Output["Output Selection"]
        Normalize{"normalize<br/>= True?"}
        NMICalc["NMI = MI / H_XY<br/>Replace NaN with 0"]
        Return["Return MI or NMI"]
    end
    
    Function --> LengthCheck
    LengthCheck -->|Pass| BinMin
    BinMin -->|Pass| BinMax
    BinMax -->|Pass| Weighted
    
    LengthCheck -->|Fail| Error1["SSException:<br/>Length mismatch"]
    BinMin -->|Fail| Error2["SSException:<br/>Bin minimum too large"]
    BinMax -->|Fail| Error3["SSException:<br/>Bin maximum too small"]
    
    Weighted --> Hist
    Hist --> Shannon
    Shannon --> MICalc
    MICalc --> Normalize
    Normalize -->|True| NMICalc
    Normalize -->|False| Return
    NMICalc --> Return
```

**Sources:** [soursop/ssmutualinformation.py:26-105]()

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `X` | `np.array` | (required) | First observable array. Must be same length as `Y`. |
| `Y` | `np.array` | (required) | Second observable array. Must be same length as `X`. |
| `bins` | `np.array` | (required) | Histogram bins. Must be uniformly spaced, monotonically increasing, and span the full data range of both `X` and `Y`. |
| `weights` | `np.array` or `False` | `False` | Optional weights for each observation. If provided, must be same length as `X` and `Y`. |
| `normalize` | `bool` | `False` | If `True`, returns normalized mutual information (NMI) scaled by joint Shannon entropy. |

#### Returns

- **float**: Mutual information value (or normalized mutual information if `normalize=True`)
- Larger values indicate greater mutual information
- MI values are bin-size dependent
- NMI values are in range [0, 1] where 0 indicates independence and 1 indicates perfect correlation

#### Implementation Details

**Histogram Generation** [soursop/ssmutualinformation.py:79-86]():
```python
if weights:
    c_XY = np.histogram2d(X,Y,bins,weights=weights)[0]
    c_X = np.histogram(X,bins,weights=weights)[0]
    c_Y = np.histogram(Y,bins,weights=weights)[0]
else:
    c_XY = np.histogram2d(X,Y,bins)[0]
    c_X = np.histogram(X,bins)[0]
    c_Y = np.histogram(Y,bins)[0]
```

**MI Calculation** [soursop/ssmutualinformation.py:88-92]():
```python
H_X = shan_entropy(c_X)
H_Y = shan_entropy(c_Y)
H_XY = shan_entropy(c_XY)
MI = H_X + H_Y - H_XY
```

**Normalization** [soursop/ssmutualinformation.py:94-103]():
```python
if normalize:
    nmi = MI / H_XY
    if np.isnan(nmi):
        nmi = 0
    return nmi
```

**Sources:** [soursop/ssmutualinformation.py:26-105]()

---

### shan_entropy()

Calculate the Shannon entropy of a probability distribution.

```mermaid
graph LR
    Input["c: np.array<br/>(1D or 2D histogram counts)"]
    
    subgraph Processing["Entropy Calculation"]
        Normalize["c_normalized = c / sum(c)"]
        NonZero["Extract non-zero elements"]
        Entropy["H = -Σ(c_i * log(c_i))"]
    end
    
    Output["Shannon Entropy H"]
    
    Input --> Normalize
    Normalize --> NonZero
    NonZero --> Entropy
    Entropy --> Output
    
    style Processing fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:110-136]()

#### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `c` | `np.array` | Array (1D or 2D) or list of numerical values representing histogram counts or probability distribution |

#### Returns

- **float**: Shannon entropy of the distribution
- Higher entropy indicates more uniform (less predictable) distribution
- Lower entropy indicates more concentrated (more predictable) distribution

#### Implementation

[soursop/ssmutualinformation.py:127-135]():
```python
# normalize such that all elements sum up to 1
c_normalized = c / float(np.sum(c))

# now convert into a single vector of non-zero elements
c_normalized = c_normalized[np.nonzero(c_normalized)]

# compute the entropy associated with this vector
H = -sum(c_normalized* np.log(c_normalized))
return H
```

**Sources:** [soursop/ssmutualinformation.py:110-136]()

---

## Usage Patterns

### Pattern 1: Correlation Between Global Properties

Analyze the relationship between radius of gyration and end-to-end distance across a trajectory.

```mermaid
sequenceDiagram
    participant User
    participant SSProtein
    participant ssmutualinformation
    
    User->>SSProtein: protein = traj.protein(0)
    User->>SSProtein: rg = protein.get_radius_of_gyration()
    User->>SSProtein: e2e = protein.get_end_to_end_distance()
    User->>User: bins = np.arange(min, max, step)
    User->>ssmutualinformation: mi = calc_MI(rg, e2e, bins)
    ssmutualinformation-->>User: Mutual information value
    
    Note over User,ssmutualinformation: Higher MI indicates Rg and E2E<br/>are strongly correlated
```

**Sources:** [soursop/ssmutualinformation.py:26-105]()

---

### Pattern 2: Weighted Analysis

Apply frame-specific weights (e.g., from reweighting schemes) to the mutual information calculation.

```mermaid
graph TB
    subgraph Data["Trajectory Data"]
        Frames["N frames"]
        Obs1["Observable X:<br/>per-frame values"]
        Obs2["Observable Y:<br/>per-frame values"]
        Weights["Frame weights:<br/>from reweighting"]
    end
    
    subgraph Analysis["Weighted MI Analysis"]
        Bins["Define bins spanning<br/>data range"]
        CalcMI["calc_MI(X, Y, bins,<br/>weights=weights)"]
    end
    
    subgraph Result["Interpretation"]
        MI["Weighted MI<br/>accounts for frame importance"]
    end
    
    Frames --> Obs1
    Frames --> Obs2
    Frames --> Weights
    
    Obs1 --> CalcMI
    Obs2 --> CalcMI
    Weights --> CalcMI
    Bins --> CalcMI
    
    CalcMI --> MI
    
    style Analysis fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:48-51](), [soursop/ssmutualinformation.py:79-86]()

---

### Pattern 3: Normalized MI for Comparing Observables

Use normalized mutual information to compare correlation strength across different observable pairs.

| Observable Pair | MI (Raw) | NMI (Normalized) | Interpretation |
|----------------|----------|------------------|----------------|
| Rg vs E2E | 0.693 | 0.95 | Very strong correlation |
| Rg vs SASA | 0.450 | 0.65 | Moderate correlation |
| E2E vs χ1 | 0.120 | 0.15 | Weak correlation |

```mermaid
graph LR
    subgraph Pairs["Observable Pairs"]
        Pair1["X1, Y1"]
        Pair2["X2, Y2"]
        Pair3["X3, Y3"]
    end
    
    subgraph Calculate["Calculate NMI"]
        NMI1["calc_MI(X1, Y1, bins,<br/>normalize=True)"]
        NMI2["calc_MI(X2, Y2, bins,<br/>normalize=True)"]
        NMI3["calc_MI(X3, Y3, bins,<br/>normalize=True)"]
    end
    
    subgraph Compare["Compare Results"]
        Scale["All values in [0, 1]<br/>regardless of bin size"]
        Rank["Rank correlations<br/>by NMI magnitude"]
    end
    
    Pair1 --> NMI1
    Pair2 --> NMI2
    Pair3 --> NMI3
    
    NMI1 --> Scale
    NMI2 --> Scale
    NMI3 --> Scale
    
    Scale --> Rank
    
    style Calculate fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:53-60](), [soursop/ssmutualinformation.py:94-103]()

---

## Practical Examples

### Example 1: Basic MI Calculation

```python
import numpy as np
from soursop import ssmutualinformation

# Extract observables from trajectory
# (assuming 'protein' is an SSProtein object)
rg_values = protein.get_radius_of_gyration()
e2e_values = protein.get_end_to_end_distance()

# Define bins that span both observables
bins = np.arange(0, 50, 0.5)  # 0 to 50 Å in 0.5 Å steps

# Calculate mutual information
mi = ssmutualinformation.calc_MI(rg_values, e2e_values, bins)
print(f"Mutual Information: {mi:.3f}")
```

**Sources:** [soursop/ssmutualinformation.py:26-68]()

---

### Example 2: Normalized MI for Scale-Independent Comparison

```python
# Compare multiple observable pairs
observables = {
    'Rg vs E2E': (rg_values, e2e_values),
    'Rg vs SASA': (rg_values, sasa_values),
    'E2E vs Asphericity': (e2e_values, asph_values)
}

results = {}
for name, (X, Y) in observables.items():
    # Use appropriate bins for each pair
    bins = np.linspace(min(X.min(), Y.min()), 
                       max(X.max(), Y.max()), 50)
    
    # Calculate normalized MI
    nmi = ssmutualinformation.calc_MI(X, Y, bins, normalize=True)
    results[name] = nmi

# Rank by correlation strength
for name, nmi in sorted(results.items(), key=lambda x: x[1], reverse=True):
    print(f"{name}: NMI = {nmi:.3f}")
```

**Sources:** [soursop/ssmutualinformation.py:94-103]()

---

### Example 3: Weighted MI with Reweighting

```python
# Apply frame weights from a reweighting scheme
frame_weights = compute_reweighting_factors(protein)  # external function

# Calculate weighted MI
mi_weighted = ssmutualinformation.calc_MI(
    rg_values, 
    e2e_values, 
    bins, 
    weights=frame_weights
)

# Compare to unweighted
mi_unweighted = ssmutualinformation.calc_MI(rg_values, e2e_values, bins)

print(f"Unweighted MI: {mi_unweighted:.3f}")
print(f"Weighted MI:   {mi_weighted:.3f}")
```

**Sources:** [soursop/ssmutualinformation.py:79-86]()

---

## Error Handling and Edge Cases

### Validation Checks

```mermaid
graph TD
    Start["calc_MI() called"]
    
    Check1{"len(X) == len(Y)?"}
    Check2{"min(bins) ≤<br/>min(X, Y)?"}
    Check3{"max(bins) ≥<br/>max(X, Y)?"}
    Check4{"NMI calculation<br/>produces NaN?"}
    
    Error1["SSException:<br/>X and Y vectors must be<br/>the same length"]
    Error2["SSException:<br/>Bin minimum exceeds<br/>data minimum"]
    Error3["SSException:<br/>Bin maximum below<br/>data maximum"]
    
    Replace["Replace NaN with 0"]
    Success["Return MI/NMI"]
    
    Start --> Check1
    Check1 -->|No| Error1
    Check1 -->|Yes| Check2
    Check2 -->|No| Error2
    Check2 -->|Yes| Check3
    Check3 -->|No| Error3
    Check3 -->|Yes| Check4
    Check4 -->|Yes| Replace
    Check4 -->|No| Success
    Replace --> Success
    
    style Error1 fill:#ffe6e6
    style Error2 fill:#ffe6e6
    style Error3 fill:#ffe6e6
```

**Sources:** [soursop/ssmutualinformation.py:70-77](), [soursop/ssmutualinformation.py:100-102]()

### Common Exceptions

| Exception | Condition | Line Reference |
|-----------|-----------|----------------|
| `SSException` | `len(X) != len(Y)` | [soursop/ssmutualinformation.py:70-71]() |
| `SSException` | `min(bins) > min(X)` or `min(bins) > min(Y)` | [soursop/ssmutualinformation.py:73-74]() |
| `SSException` | `max(bins) < max(X)` or `max(bins) < max(Y)` | [soursop/ssmutualinformation.py:76-77]() |

**NaN Handling:** When normalization produces NaN (typically when joint entropy is zero), the value is automatically replaced with 0 [soursop/ssmutualinformation.py:100-102]().

**Sources:** [soursop/ssmutualinformation.py:70-77](), [soursop/ssmutualinformation.py:100-102]()

---

## Testing and Validation

The module includes comprehensive tests demonstrating expected behavior.

### Test Cases

```mermaid
graph TB
    subgraph Tests["Test Suite"]
        Test1["test_shan_entropy():<br/>Validate entropy calculation<br/>for uniform distributions"]
        Test2["test_calc_MI():<br/>Perfect correlation:<br/>calc_MI(X, X, bins) ≈ 0.693"]
        Test3["test_calc_NMI():<br/>Independence: NMI(X, Y) ≈ 0<br/>Perfect: NMI(X, X) ≈ 1.0"]
    end
    
    subgraph Expected["Expected Results"]
        E1["Uniform dist. of size 2:<br/>H = 0.693"]
        E2["Uniform dist. of size 3:<br/>H = 1.099"]
        E3["Uniform dist. of size 4:<br/>H = 1.386"]
        E4["Independent vars:<br/>NMI ≈ 0"]
        E5["Identical vars:<br/>NMI ≈ 1.0"]
    end
    
    Test1 --> E1
    Test1 --> E2
    Test1 --> E3
    Test2 --> E5
    Test3 --> E4
    Test3 --> E5
    
    style Tests fill:#f9f9f9
```

**Sources:** [soursop/tests/test_ssmutual_information.py:1-31]()

### Test Implementation

**Shannon Entropy Tests** [soursop/tests/test_ssmutual_information.py:6-11]():
- Validates entropy for uniform distributions of different sizes
- Verifies logarithmic growth with distribution size

**Mutual Information Tests** [soursop/tests/test_ssmutual_information.py:14-21]():
- Tests perfect correlation: MI between a variable and itself
- Uses random normal distributions to ensure robustness

**Normalized MI Tests** [soursop/tests/test_ssmutual_information.py:22-29]():
- Tests independence: random independent variables should have NMI ≈ 0
- Tests perfect correlation: identical variables should have NMI ≈ 1.0

**Sources:** [soursop/tests/test_ssmutual_information.py:1-31]()

---

## Implementation Notes

### Dependencies

```mermaid
graph LR
    subgraph External["External Libraries"]
        NumPy["numpy:<br/>Array operations,<br/>histogram generation"]
    end
    
    subgraph Internal["SOURSOP Modules"]
        SSExc["ssexceptions:<br/>SSWarning, SSException"]
    end
    
    subgraph Module["ssmutualinformation.py"]
        Functions["calc_MI()<br/>shan_entropy()"]
    end
    
    NumPy --> Functions
    SSExc --> Functions
    
    style Module fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:20-21]()

### Design Rationale

The module is implemented as **standalone functions** rather than a class because:

1. **Generality**: MI calculations are pure statistical operations that don't require state
2. **Flexibility**: Functions can be applied to any pair of observables from any source
3. **Simplicity**: No need to instantiate objects for one-off calculations
4. **Reusability**: Functions can be easily imported and used in different contexts

**Sources:** [soursop/ssmutualinformation.py:14-18]()

---

## Performance Considerations

### Bin Size Selection

Mutual information values are **bin-size dependent**. Consider these guidelines:

| Factor | Recommendation |
|--------|---------------|
| **Number of bins** | 30-100 bins typically provide good resolution |
| **Bin width** | Should be smaller than meaningful feature scales |
| **Data coverage** | Bins must span full data range (enforced by validation) |
| **Sample size** | More frames allow finer binning |

### Computational Cost

```mermaid
graph LR
    subgraph Operations["Computational Steps"]
        Hist["Histogram generation:<br/>O(N) per histogram"]
        Entropy["Entropy calculation:<br/>O(bins) per entropy"]
        Total["Total: O(N + bins)"]
    end
    
    subgraph Scaling["Memory and Time"]
        Memory["Memory: O(N + bins²)"]
        Time["Time: Dominated by<br/>numpy histogram operations"]
    end
    
    Hist --> Total
    Entropy --> Total
    Total --> Memory
    Total --> Time
    
    style Operations fill:#f9f9f9
```

- **Time complexity**: O(N) where N is the number of frames
- **Memory complexity**: O(N + bins²) for 2D histogram
- **Bottleneck**: `np.histogram2d()` for large datasets

**Sources:** [soursop/ssmutualinformation.py:79-92]()

---

## Integration with SSProtein

While `ssmutualinformation` functions are standalone, they integrate naturally with `SSProtein` analysis methods:

```mermaid
graph TB
    subgraph SSProtein_Methods["SSProtein Analysis Methods"]
        Global["Global Properties:<br/>get_radius_of_gyration()<br/>get_end_to_end_distance()<br/>get_asphericity()"]
        
        Local["Local Properties:<br/>get_inter_residue_distances()<br/>get_dihedral_angles()<br/>get_SASA()"]
        
        SecStruct["Secondary Structure:<br/>get_secondary_structure_DSSP()<br/>get_secondary_structure_BBSEG()"]
    end
    
    subgraph Observables["Observable Extraction"]
        Array1["Observable X<br/>(np.array)"]
        Array2["Observable Y<br/>(np.array)"]
    end
    
    subgraph MI_Analysis["ssmutualinformation"]
        CalcMI["calc_MI(X, Y, bins)"]
        Result["MI or NMI value"]
    end
    
    Global --> Array1
    Global --> Array2
    Local --> Array1
    Local --> Array2
    SecStruct --> Array1
    SecStruct --> Array2
    
    Array1 --> CalcMI
    Array2 --> CalcMI
    CalcMI --> Result
    
    style MI_Analysis fill:#f9f9f9
```

**Sources:** [soursop/ssmutualinformation.py:14-18]()

---