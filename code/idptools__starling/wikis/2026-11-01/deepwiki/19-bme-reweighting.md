# BME Reweighting

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [.gitignore](.gitignore)
- [demos/bme_reweighting_example.ipynb](demos/bme_reweighting_example.ipynb)
- [demos/theta_scan_rg.pdf](demos/theta_scan_rg.pdf)
- [demos/theta_scan_rg.png](demos/theta_scan_rg.png)
- [docs/usage/possible_issues.rst](docs/usage/possible_issues.rst)
- [starling/inference/__init__.py](starling/inference/__init__.py)
- [starling/structure/bme.py](starling/structure/bme.py)
- [starling/structure/bme_utils.py](starling/structure/bme_utils.py)
- [starling/structure/ensemble.py](starling/structure/ensemble.py)

</details>



Bayesian Maximum Entropy (BME) reweighting is a statistical framework used to refine structural ensembles by integrating experimental data. In STARLING, this is implemented as a post-processing step that adjusts the statistical weights of generated conformations to match experimental observables (e.g., Radius of Gyration, FRET, NMR) while minimizing the information-theoretic bias introduced to the original ensemble [starling/structure/bme.py:1-10]().

The implementation allows for equality constraints, as well as upper and lower bounds, providing a principled way to balance experimental fitting with the preservation of ensemble diversity [starling/structure/bme.py:58-74]().

## System Overview and Data Flow

The BME workflow typically starts with a generated `Ensemble` object. Observables are calculated for every frame in the ensemble, and these are compared against `ExperimentalObservable` targets to produce a `BMEResult` containing optimized weights.

### Conceptual to Code Mapping

| Concept | Code Entity | Role |
|:---|:---|:---|
| **Target Data** | `ExperimentalObservable` | Stores experimental value, uncertainty, and constraint type [starling/structure/bme_utils.py:24-47](). |
| **Optimization Engine** | `BME` | Performs the L-BFGS-B optimization of Lagrange multipliers [starling/structure/bme.py:108-130](). |
| **Reweighting Logic** | `Ensemble.reweight_bme()` | High-level API to perform reweighting and cache results [starling/structure/ensemble.py:488-540](). |
| **Selection Strategy** | `theta_scan` | Automates selection of the $\theta$ regularization parameter via L-curve analysis [starling/structure/bme_utils.py:270-350](). |
| **Output Container** | `BMEResult` | Holds optimized weights, final $\chi^2$, and diagnostic metrics [starling/structure/bme_utils.py:91-110](). |

**Sources:** [starling/structure/bme.py:12-22](), [starling/structure/ensemble.py:110-112]().

### BME Data Flow Diagram

```mermaid
graph TD
    subgraph "Natural Language Space"
        ExpData["Experimental Measurement (e.g. Rg = 25Å)"]
        EnsembleData["Generated Protein Structures"]
        ThetaSelection["Find balance between Fit and Diversity"]
    end

    subgraph "Code Entity Space"
        EO["ExperimentalObservable class"]
        ENS["Ensemble class"]
        CALC["Calculated values (numpy array)"]
        BME_OBJ["BME class instance"]
        SCAN["theta_scan() function"]
        RES["BMEResult dataclass"]
    end

    ExpData --> EO
    EnsembleData --> ENS
    ENS -- "ensemble.radius_of_gyration()" --> CALC
    EO & CALC --> BME_OBJ
    BME_OBJ -- "fit(theta)" --> RES
    BME_OBJ -- "L-curve analysis" --> SCAN
    SCAN --> RES
    RES -- "use_bme_weights=True" --> ENS
```
**Sources:** [starling/structure/bme.py:31-48](), [starling/structure/bme_utils.py:270-285]().

---

## Key Components

### ExperimentalObservable
This dataclass encapsulates a single experimental constraint. It supports three types of constraints defined in `VALID_CONSTRAINTS`: `"equality"`, `"upper"`, and `"lower"` [starling/structure/bme_utils.py:18-18]().

- **Equality**: The ensemble average must match the value within the specified uncertainty.
- **Upper/Lower**: Enforces bounds on the ensemble average. The implementation uses specific Lagrange multiplier bounds to enforce these: `(0.0, None)` for upper and `(None, 0.0)` for lower [starling/structure/bme_utils.py:82-85]().

**Sources:** [starling/structure/bme_utils.py:24-47]().

### The BME Class
The `BME` class handles the core optimization logic. It minimizes a global objective function:
$$\mathcal{L}(\lambda) = \theta \Gamma(\lambda) + \chi^2(\lambda)$$
where $\Gamma$ is related to the partition function and $\theta$ is the regularization parameter [starling/structure/bme.py:214-230]().

- **Initialization**: Sets up initial weights (uniform by default) and randomizes Lagrange multipliers ($\lambda$) [starling/structure/bme.py:143-151]().
- **fit()**: Uses `scipy.optimize.minimize` with the `L-BFGS-B` method to find optimal $\lambda$ values [starling/structure/bme.py:290-310]().
- **phi ($\phi$)**: A key metric representing the "effective fraction" of the ensemble remaining after reweighting, calculated as $\exp(-D_{KL})$ [starling/structure/bme_utils.py:126-133]().

**Sources:** [starling/structure/bme.py:108-160](), [starling/structure/bme_utils.py:11-15]().

### Theta ($\theta$) Regularization and L-curve
The parameter `theta` controls the trade-off between fitting experimental data (low $\theta$) and staying close to the original ensemble (high $\theta$). 

The `theta_scan` function automates the selection of an optimal $\theta$ by:
1. Performing reweighting across a range of $\theta$ values [starling/structure/bme_utils.py:290-310]().
2. Plotting the "L-curve": $\chi^2$ vs. $\phi$ (or $D_{KL}$) [starling/structure/bme_utils.py:350-380]().
3. Identifying the "elbow" of the curve where further fitting leads to a drastic loss in ensemble diversity [starling/structure/bme_utils.py:315-325]().

**Sources:** [starling/structure/bme_utils.py:270-350]().

---

## Integration with Ensemble

The `Ensemble` class provides the primary interface for BME via the `reweight_bme` method.

### Workflow Diagram

```mermaid
sequenceDiagram
    participant U as User
    participant E as Ensemble
    participant B as BME
    participant R as BMEResult

    U->>E: ensemble.reweight_bme(observables, calculated_values)
    E->>B: BME(observables, calculated_values)
    B->>B: fit(theta)
    B-->>R: Create BMEResult
    R-->>E: Return Result
    E->>E: Cache result in self.__bme_result
    U->>E: ensemble.radius_of_gyration(use_bme_weights=True)
    E->>E: Apply weights from cached BMEResult
    E-->>U: Return reweighted mean
```
**Sources:** [starling/structure/ensemble.py:488-540](), [starling/structure/ensemble.py:110-112]().

### Diagnostics and Quality Control
The `BMEResult.diagnostics()` method provides critical feedback on the reweighting quality [starling/structure/bme_utils.py:137-180]():
- **neff_entropy**: The effective sample size based on Shannon entropy ($N \cdot \phi$).
- **neff_renyi2**: A more conservative effective sample size based on the participation ratio ($1 / \sum w_i^2$).
- **chi2_improvement**: The percentage reduction in $\chi^2$ error.
- **Warnings**: Triggered if $\phi$ drops below a threshold (default 0.5), indicating potential overfitting or a poor initial ensemble [starling/structure/bme_utils.py:200-210]().

**Sources:** [starling/structure/bme_utils.py:137-210]().

---

## Implementation Details

### Objective Function
The optimization minimizes the log-partition function like objective:
```python
# Simplified logic from BME._objective_function
def objective(lambdas):
    # Calculate weights for current lambdas
    log_w = log_initial_w - (calculated_values @ lambdas) / theta
    log_z = logsumexp(log_w)
    
    # Calculate chi-squared
    current_means = weights @ calculated_values
    chi2 = sum(((current_means - exp_values) / uncertainties)**2)
    
    return theta * log_z + 0.5 * chi2
```
**Sources:** [starling/structure/bme.py:214-265]().

### Constraints Handling
During `fit()`, the `L-BFGS-B` optimizer respects the bounds provided by `ExperimentalObservable.get_bounds()`. This ensures that for an "upper" constraint, the Lagrange multiplier remains positive, only allowing the optimizer to penalize (push down) values that exceed the target [starling/structure/bme_utils.py:71-86]().

**Sources:** [starling/structure/bme.py:300-305]().

---