# Model Utilities: Precision, Tiling, Strategy

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/bioemu2/biomoa/utils/__init__.py](src/bioemu2/biomoa/utils/__init__.py)
- [src/bioemu2/biomoa/utils/feature_converter.py](src/bioemu2/biomoa/utils/feature_converter.py)
- [src/bioemu2/biomoa/utils/precision.py](src/bioemu2/biomoa/utils/precision.py)
- [src/bioemu2/biomoa/utils/strategy.py](src/bioemu2/biomoa/utils/strategy.py)
- [src/bioemu2/biomoa/utils/tiling.py](src/bioemu2/biomoa/utils/tiling.py)
- [src/bioemu2/biomoa/utils/warnings.py](src/bioemu2/biomoa/utils/warnings.py)

</details>



The `src/bioemu2/biomoa/utils/` package provides low-level supporting infrastructure for the BioMoA score model stack. These utilities handle feature transformations between dense structure representations and model-consumable flat dictionaries, manage fine-grained matrix multiplication precision knobs, instantiate model sub-modules dynamically via a centralized strategy factory, resolve memory-saving tiling and chunking configurations with automatic variant switching, and suppress expected third-party deprecation warnings.

---

## 1. Feature Conversion (`feature_converter.py`)

The BioMoA score model consumes flat feature dictionaries rather than nested object hierarchies. The module `src/bioemu2/biomoa/utils/feature_converter.py` bridges dense structure batches (`BatchAllAtomStructure`) to flat tensor dictionaries expected by the transformer trunk and its sub-modules [src/bioemu2/biomoa/utils/feature_converter.py:1-24].

### Key Functions and Behavior
- **`convert_densebatch_to_dictionary(dense_batch, *, build_pair_features=True)`**: Extracts zero-copy tensor references from `dense_batch.to_dict()` and transforms specific keys (boolean casts, derived token features, padding sentinels, and contact conditioning fallbacks) [src/bioemu2/biomoa/utils/feature_converter.py:19-51].
- **Index Sentinels**: Uses `INDEX_PAD_VALUE: int = -1` for index-valued fields to ensure unmasked padding values fail loudly via scatter out-of-bounds errors rather than silently aliasing index `0` [src/bioemu2/biomoa/utils/feature_converter.py:13-16, 86].
- **Contact Conditioning Fallback**: If `build_pair_features` is active and contact conditioning is absent, it invokes `make_unspecified_contact_conditioning` and `make_unspecified_contact_threshold` [src/bioemu2/biomoa/utils/feature_converter.py:58-66].
- **Derived Features**: Gathers residue types and indices from atom-level tensors via `token_center_atom_index`, computes symmetric `token_bonds`, and transposes MSA features to depth-first layout [src/bioemu2/biomoa/utils/feature_converter.py:75-120].

```mermaid
graph TD
    A["BatchAllAtomStructure"] --> B["convert_densebatch_to_dictionary"]
    B --> C["TensorReferenceExtraction"]
    B --> D["BooleanCastMasks"]
    B --> E["ContactConditioningFallback"]
    B --> F["DerivedTokenFeatures"]
    B --> G["MSATransposeAndPadding"]
    C --> H["FlatFeatureDictionary"]
    D --> H
    E --> H
    F --> H
    G --> H

    sub-sources-1
    Sources: src/bioemu2/biomoa/utils/feature_converter.py:1-122
```

Sources: `src/bioemu2/biomoa/utils/feature_converter.py:1-122`

---

## 2. Precision Configuration and Matmul Knobs (`precision.py`)

Matrix multiplication precision in kernel sub-modules (such as outer product mean and triangle multiplication) is governed by `src/bioemu2/biomoa/utils/precision.py` [src/bioemu2/biomoa/utils/precision.py:1-19].

### Precision Rules and Resolution
- **Knob Values**: Supports `"ieee"`, `"tf32"`, and `"tf32x3"` (with `cueq` TriMul also permitting `None`) [src/bioemu2/biomoa/utils/precision.py:1-5].
- **Variant Validation**: Maps specific precision knobs to compatible execution variants via `VARIANTS_WITH_PRECISION_KNOB` and `DEFAULT_PRECISION_KEYS` [src/bioemu2/biomoa/utils/precision.py:26-41].
- **`resolve_precision_cfg(precision_cfg, *, site, variant)`**: Traverses a nested precision configuration dictionary (e.g., `msa_module.opm_matmul_precision`), falling back to variant defaults if `null` is specified, and warns if an incompatible variant attempts to claim a knob without switching variants [src/bioemu2/biomoa/utils/precision.py:44-85].

Sources: `src/bioemu2/biomoa/utils/precision.py:1-85`

---

## 3. Strategy Factory and Dynamic Instantiation (`strategy.py`)

The strategy module `src/bioemu2/biomoa/utils/strategy.py` supplies a centralized factory for instantiating variant-based sub-modules [src/bioemu2/biomoa/utils/strategy.py:1-5].

### The `make()` Factory
- **`make(registry, variant, *, matmul_precision=None, **kwargs)`**: Looks up `variant` (string or enum) in a module registry dictionary, dynamically imports the target module path, inspects the class signature for `matmul_precision`, injects it if supported, and passes remaining `**kwargs` to the constructor [src/bioemu2/biomoa/utils/strategy.py:13-50].
- **Silent Ignoring**: Variants that do not accept a precision knob (e.g., vanilla or transition modules) silently ignore `matmul_precision`, enabling uniform parameter passing across kernel families [src/bioemu2/biomoa/utils/strategy.py:21-27].

```mermaid
graph TD
    A["CallerConfig"] -->|RegistryAndVariant| B["strategy.make"]
    B --> C["LookupRegistryEntry"]
    C --> D["ImportModuleAndClass"]
    D --> E["InspectConstructorSignature"]
    E -->|SupportsMatmulPrecision| F["InjectMatmulPrecision"]
    E -->|DoesNotSupport| G["SkipPrecisionKwarg"]
    F --> H["InstantiateModule"]
    G --> H
    H --> I["nn.ModuleInstance"]

    sub-sources-3
    Sources: src/bioemu2/biomoa/utils/strategy.py:1-50
```

Sources: `src/bioemu2/biomoa/utils/strategy.py:1-50`

---

## 4. Memory Tiling and Chunking (`tiling.py`)

To control memory consumption during long-sequence or deep-MSA inference, `src/bioemu2/biomoa/utils/tiling.py` parses, validates, and auto-corrects tiling configurations (`tiling_cfg`) [src/bioemu2/biomoa/utils/tiling.py:1-32].

### Key Data Structures and Functions
- **Module Known Keys**: `MODULE_KNOWN_KEYS` defines which tiling keys are legal for specific consumers (e.g., `BioMoaCModel` vs `PairMixerSequenceEmbedder`) [src/bioemu2/biomoa/utils/tiling.py:61-64].
- **Global Dominance**: Global keys (`tile_size` and `chunk_size`) override all corresponding `*_tile_size` or `*_chunk_size` entries [src/bioemu2/biomoa/utils/tiling.py:39-42, 131-165].
- **Variant Auto-Switching**: `resolve_tiling_cfg` checks `TILING_KEY_VARIANT_REQUIREMENTS`. If a non-default tile or chunk size is requested on an incompatible variant, it automatically switches the variant to the corresponding entry in `PREFERRED_VARIANT_FOR_TILING` (e.g., switching to `triton_memory_saving` for `opm_tile_size`) and emits a warning [src/bioemu2/biomoa/utils/tiling.py:69-94, 131-171].

Sources: `src/bioemu2/biomoa/utils/tiling.py:1-171`

---

## 5. Warning Management (`warnings.py`)

Third-party dependencies (such as `cuequivariance`) occasionally trigger deprecation warnings at import time that cannot be directly modified within BioEmu2. The module `src/bioemu2/biomoa/utils/warnings.py` isolates these warnings [src/bioemu2/biomoa/utils/warnings.py:1-28].

- **`suppress_cuequivariance_jit_deprecation()`**: A context manager that catches `torch.jit.script` deprecation warnings matching `_TORCH_JIT_SCRIPT_DEPRECATIONS` and suppresses them during third-party library initialization [src/bioemu2/biomoa/utils/warnings.py:5-28].

Sources: `src/bioemu2/biomoa/utils/warnings.py:1-28`

---