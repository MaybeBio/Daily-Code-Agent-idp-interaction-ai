# Search Engine

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [starling/configs.py](starling/configs.py)
- [starling/inference/benchmark_mds.py](starling/inference/benchmark_mds.py)
- [starling/scripts/starling_pretokenize.py](starling/scripts/starling_pretokenize.py)
- [starling/scripts/starling_search.py](starling/scripts/starling_search.py)
- [starling/search/__init__.py](starling/search/__init__.py)
- [starling/search/builder.py](starling/search/builder.py)
- [starling/search/search_engine.py](starling/search/search_engine.py)
- [starling/search/search_utils.py](starling/search/search_utils.py)
- [starling/search/similarity_search.py](starling/search/similarity_search.py)
- [starling/search/store.py](starling/search/store.py)

</details>



The STARLING Search Engine provides a high-performance similarity search system for protein ensembles. Built on top of **FAISS** (Facebook AI Similarity Search), it enables billion-scale Approximate Nearest Neighbor (ANN) lookups, multi-level metadata filtering, and exact sequence reranking. The system bridges the gap between latent embedding space and raw sequence data by coupling FAISS indexes with a compressed SQLite-backed metadata store.

### System Overview

The search infrastructure is composed of three primary software layers:

1.  **Storage Layer**: Consists of the serialized FAISS index (e.g., `.faiss`) and the `SequenceStore` (SQLite), which holds raw sequences, headers, and pre-computed hashes [starling/search/store.py:22-34]().
2.  **Engine Layer**: The `SearchEngine` class coordinates ANN lookups, length-based pre-filtering, and a chain of `CandidateFilter` objects [starling/search/search_engine.py:94-100]().
3.  **Tooling Layer**: CLI tools like `starling-search` and `starling-pretokenize` for index construction and querying [starling/scripts/starling_search.py:27-33]().

### Component Interaction

The following diagram illustrates the flow from raw sequence data to a searchable index and the subsequent query process.

**Search System Data Flow**
```mermaid
graph TD
    subgraph "Build Phase"
        FASTA["FASTA Files"] --> SPT["starling-pretokenize"]
        SPT --> TOKENS["Tokenized PT Files"]
        FEAT["Feature Shards (.pt)"] --> IB["IndexBuilder.build_index()"]
        TOKENS --> IB
        IB --> FAISS_IDX["FAISS Index (.faiss)"]
        IB --> SEQ_DB["SequenceStore (.seqs.sqlite)"]
    end

    subgraph "Query Phase"
        QUERY_SEQ["Query Sequence"] --> SE_ENC["sequence_encoder_backend"]
        SE_ENC --> EMB["Query Embedding"]
        EMB --> SE_SEARCH["SearchEngine.search()"]
        FAISS_IDX --> SE_SEARCH
        SEQ_DB --> SE_SEARCH
        SE_SEARCH --> FILTERS["CandidateFilter Chain"]
        FILTERS --> RESULTS["Ranked Candidates"]
    end
```
**Sources:** [starling/search/builder.py:22-39](), [starling/search/search_engine.py:9-15](), [starling/scripts/starling_pretokenize.py:2-11]().

---

### Key Entities and Architecture

The search system relies on specific classes to manage the lifecycle of an index.

| Entity | File Path | Role |
| :--- | :--- | :--- |
| `SearchEngine` | [starling/search/search_engine.py:93]() | Main interface for executing queries and applying filters. |
| `IndexBuilder` | [starling/search/builder.py:152]() | Orchestrates FAISS training (OPQ+IVF-PQ) and shard discovery. |
| `SequenceStore` | [starling/search/store.py:215]() | Manages SQLite storage for sequences with zstd compression. |
| `Candidate` | [starling/search/search_utils.py:220]() | Data container for a single search hit (score, gid, meta). |
| `CandidateFilter` | [starling/search/search_utils.py:212]() | Abstract base class for post-ANN filtering logic. |

**Code Entity Space Mapping**
```mermaid
graph LR
    subgraph "starling.search"
        SE["SearchEngine"] -- "uses" --> SS["SequenceStore"]
        SE -- "uses" --> SC["ScoreConverter"]
        IB["IndexBuilder"] -- "creates" --> SS
        IB -- "trains" --> FAISS["faiss.Index"]
    end

    subgraph "starling.search.search_utils"
        SC -- "produces" --> CAND["Candidate"]
        CAND -- "passed to" --> CF["CandidateFilter"]
        CF --> LF["LengthFilter"]
        CF --> EF["ExactMatchFilter"]
        CF --> IF["SequenceIdentityFilter"]
    end
```
**Sources:** [starling/search/search_engine.py:72-83](), [starling/search/builder.py:124-127](), [starling/search/search_utils.py:72-81]().

---

### Core Functionality

#### Index Construction and SequenceStore
The construction process involves training a FAISS index using Optimized Product Quantization (OPQ) and Inverted File (IVF) clusters. This allows the system to scale to millions of protein sequences while maintaining a small memory footprint by compressing vectors to Product Quantization (PQ) codes [starling/search/builder.py:10-16](). Simultaneously, a `SequenceStore` is built to provide O(log N) access to metadata like sequence length and 8-byte hashes for deduplication [starling/search/store.py:9-16]().

For details, see [Index Construction and SequenceStore](#8.1).

#### Querying and Filtering
Querying is a two-stage process. First, an ANN search retrieves coarse candidates from the FAISS index. Second, these candidates are passed through a chain of filters (e.g., `LengthFilter`, `ExactMatchFilter`) to ensure they meet user-specified constraints [starling/search/search_engine.py:9-13](). The system also supports reranking, where top candidates are re-scored using the full encoder for higher precision [starling/search/search_engine.py:38-39]().

For details, see [Querying and Filtering](#8.2).

---

### CLI Tools

The search engine is primarily accessed via two command-line utilities:

*   **`starling-pretokenize`**: Processes FASTA files into tokenized `.pt` files required for building the sequence store [starling/scripts/starling_pretokenize.py:2-6]().
*   **`starling-search`**: A unified tool for building new indexes (`build` command) and querying existing ones (`query` command) [starling/scripts/starling_search.py:27-33]().

**Sources:** [starling/scripts/starling_pretokenize.py:44-55](), [starling/scripts/starling_search.py:47-105]().

---