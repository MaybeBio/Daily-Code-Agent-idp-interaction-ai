# Input Format and Preparation

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [example/output.tar.gz](example/output.tar.gz)
- [example/pae_model_4_multimer_v3_pred_4.json](example/pae_model_4_multimer_v3_pred_4.json)
- [example/test.csv](example/test.csv)
- [example/test.fasta](example/test.fasta)
- [example/unrelaxed_model_4_multimer_v3_pred_4.pdb](example/unrelaxed_model_4_multimer_v3_pred_4.pdb)
- [run_disobind.py](run_disobind.py)

</details>



This page describes how to prepare input data for Disobind predictions, including the CSV input format, UniProt ID requirements, sequence download procedures, and embedding generation. For information about running predictions and model configuration, see [Running Predictions](). For details on AlphaFold integration with predictions, see [AlphaFold Integration]().

## Overview

Disobind accepts protein pairs as input via a CSV file. Each row specifies one protein fragment pair for prediction. The system automatically downloads UniProt sequences, generates T5 embeddings, and prepares input tensors for the neural network. Input preparation is handled by the `Disobind` class in `run_disobind.py` [run_disobind.py:44-106]() and the `Embeddings` class in `dataset/create_input_embeddings.py` [dataset/create_input_embeddings.py:38-40]().

**Sources:** [README.md:42-64](), [run_disobind.py:44-106]()

## Input CSV Format

### Basic Format: Disobind Only

For standard Disobind predictions without AlphaFold integration, each CSV row contains 6 comma-separated fields [run_disobind.py:7-8]():

```
UniProt_ID1,start1,end1,UniProt_ID2,start2,end2
```

**Field Descriptions:**

| Field | Type | Description |
|-------|------|-------------|
| `UniProt_ID1` | String | UniProt accession for protein 1 (must be an IDR) |
| `start1` | Integer | Starting residue position in UniProt sequence for protein 1 |
| `end1` | Integer | Ending residue position in UniProt sequence for protein 1 |
| `UniProt_ID2` | String | UniProt accession for protein 2 (may or may not be an IDR) |
| `start2` | Integer | Starting residue position in UniProt sequence for protein 2 |
| `end2` | Integer | Ending residue position in UniProt sequence for protein 2 |

**Example:**
```csv
P04273,95,193,P04273,95,192
```

**Sources:** [README.md:48-58](), [run_disobind.py:208-255](), [example/test.csv:1]()

### Extended Format: Disobind + AlphaFold2/3

For combining Disobind predictions with AlphaFold structural predictions, each CSV row contains 12 comma-separated fields [run_disobind.py:9-10]():

```
UniProt_ID1,start1,end1,UniProt_ID2,start2,end2,AF2_struct_file_path,AF2_pae_file_path,chain1,chain2,offset1,offset2
```

**Additional Field Descriptions:**

| Field | Type | Description |
|-------|------|-------------|
| `AF2_struct_file_path` | String | Path to AlphaFold structure file (.pdb or .cif) |
| `AF2_pae_file_path` | String | Path to AlphaFold PAE JSON file |
| `chain1` | String | Chain ID in AF2 structure corresponding to protein 1 fragment |
| `chain2` | String | Chain ID in AF2 structure corresponding to protein 2 fragment |
| `offset1` | Integer | Residue position offset between AF2 structure and UniProt numbering for protein 1 |
| `offset2` | Integer | Residue position offset between AF2 structure and UniProt numbering for protein 2 |

**Offset Calculation:** The offset is the difference between the AF2 structure residue numbering and the UniProt residue numbering [README.md:63-65](). Set to 0 if the AF2 structure corresponds to the full UniProt sequence or the exact fragment specified.

**Example:**
```csv
P04273,95,193,P04273,95,193,./example/unrelaxed_model_4_multimer_v3_pred_4.pdb,./example/pae_model_4_multimer_v3_pred_4.json,B,C,0,0
```

**Sources:** [README.md:57-65](), [run_disobind.py:232-253](), [example/test.csv:2]()

## Input Requirements and Constraints

The following requirements and constraints apply to all input:

1. **Binary Complexes Only:** Disobind only processes binary (two-protein) complexes. For non-binary complexes ABC, create binary pairs: AB, BC, AC [README.md:43-44]().
2. **Protein 1 Must Be an IDR:** The first protein in each pair must be an intrinsically disordered region. Protein 2 may be ordered or disordered [README.md:45]().
3. **Interacting Pairs:** The input protein pair is assumed to be interacting. Disobind predicts *where* they interact, not *whether* they interact [README.md:44]().
4. **Valid UniProt IDs:** By default, CSV input expects proteins with UniProt accessions for sequence retrieval [README.md:49]().
5. **Alternative FASTA Input:** For proteins without UniProt accessions, the FASTA format can be used to provide sequences directly [README.md:50-51]().

**Sources:** [README.md:42-51](), [run_disobind.py:4-10]()

## Input Processing Pipeline

The following diagram shows the complete input processing workflow from CSV file to model-ready tensors:

**Diagram: Input Processing Workflow**

```mermaid
flowchart TD
    CSVInput["CSV Input File"]
    ReadCSV["process_input_file()<br/>[run_disobind.py:212]"]
    ParseRows["Parse CSV Rows<br/>Extract UniProt IDs<br/>and Residue Ranges"]
    EntryDict["entry_ids List<br/>af_dict Dictionary"]
    
    CheckUniSeq{"UniProt Sequences<br/>Cached?"}
    UniSeqFile["UniProt_seq.json<br/>[output_dir/]"]
    DownloadSeq["download_uniprot_seq()<br/>[run_disobind.py:299]"]
    GetUniqueIDs["get_unique_uni_ids()<br/>[run_disobind.py:280]"]
    ParallelDownload["parallelize_uni_seq_download()<br/>[run_disobind.py:261]<br/>Uses multiprocessing.Pool"]
    APICall["get_uniprot_seq()<br/>[from_APIs_with_love.py:10]"]
    
    ProcessPairs["process_input_pairs()<br/>[run_disobind.py:333]"]
    ValidatePairs["Validate UniProt IDs<br/>in Downloaded Sequences"]
    ProtPairs["prot_pairs List<br/>Valid Entry IDs"]
    
    CreateFasta["create_fasta_from_headers()<br/>[create_input_embeddings.py:186]"]
    FastaFile["FASTA File<br/>p1_p2_test.fasta"]
    
    CheckEmb{"Embeddings<br/>Exist?"}
    EmbFile["Embeddings HDF5<br/>p1_p2_test.h5"]
    GenEmb["get_embeddings()<br/>[utility.py:155]"]
    T5Script["ProtT5_embeddings()<br/>[utility.py:342]"]
    
    LoadEmb["Load Embeddings<br/>get_global_embeddings()<br/>[create_input_embeddings.py:276]"]
    EmbDict["prot1_emb, prot2_emb<br/>Dictionaries"]
    
    BatchProcess["get_predictions()<br/>[run_disobind.py:168]<br/>Process in Batches"]
    PrepTensors["get_input_tensors()<br/>[run_disobind.py:475]"]
    Padding["apply_padding()<br/>Pad to 100x100"]
    CoarseGrain["prepare_input()<br/>[src/utils.py:10]<br/>Apply Coarse-Graining"]
    
    ModelInput["Model-Ready Tensors<br/>prot1, prot2, target_mask"]
    
    CSVInput --> ReadCSV
    ReadCSV --> ParseRows
    ParseRows --> EntryDict
    
    EntryDict --> CheckUniSeq
    CheckUniSeq -->|No| DownloadSeq
    CheckUniSeq -->|Yes| ProcessPairs
    
    DownloadSeq --> GetUniqueIDs
    GetUniqueIDs --> ParallelDownload
    ParallelDownload --> APICall
    APICall --> UniSeqFile
    UniSeqFile --> ProcessPairs
    
    ProcessPairs --> ValidatePairs
    ValidatePairs --> ProtPairs
    ProtPairs --> CreateFasta
    CreateFasta --> FastaFile
    
    FastaFile --> CheckEmb
    CheckEmb -->|No| GenEmb
    CheckEmb -->|Yes| LoadEmb
    
    GenEmb --> T5Script
    T5Script --> EmbFile
    EmbFile --> LoadEmb
    
    LoadEmb --> EmbDict
    EmbDict --> BatchProcess
    BatchProcess --> PrepTensors
    PrepTensors --> Padding
    Padding --> CoarseGrain
    CoarseGrain --> ModelInput
```

**Sources:** [run_disobind.py:111-207](), [run_disobind.py:212-368](), [run_disobind.py:475-529](), [dataset/create_input_embeddings.py:186-273]()

## CSV Parsing and Entry ID Format

The `process_input_file()` method (specifically for CSV type) parses the file and creates internal entry identifiers [run_disobind.py:212-255]():

```
UniID1:start1:end1--UniID2:start2:end2_0
```

The `_0` suffix indicates this is the first fragment pair for this combination. The method also constructs an `af_dict` dictionary to store AlphaFold-related information.

**Diagram: CSV Parsing and Data Structures**

```mermaid
flowchart LR
    subgraph "CSV Input"
        Row["P04273,95,193,P04273,95,192<br/>OR<br/>P04273,95,193,P04273,95,193,<br/>path/to/pdb,path/to/json,<br/>B,C,0,0"]
    end
    
    subgraph "read_csv() Processing"
        Split["Split by Comma<br/>Validate Length:<br/>6 fields or 12 fields"]
        ExtractFields["Extract Fields:<br/>uni_id1, start1, end1<br/>uni_id2, start2, end2<br/>(+AF2 fields if present)"]
        CreateID["Create Entry ID:<br/>entry_id = f'{uni_id1}:{start1}:{end1}--<br/>{uni_id2}:{start2}:{end2}_0'"]
    end
    
    subgraph "Data Structures"
        EntryList["entry_ids List<br/>['P04273:95:193--P04273:95:192_0']"]
        AFDict["af_dict Dictionary<br/>{<br/>  entry_id: {<br/>    'struct_file': path,<br/>    'json_file': path,<br/>    'required_chains': {<br/>      'chains': [chain1, chain2],<br/>      'offsets': [offset1, offset2]<br/>    }<br/>  }<br/>}"]
    end
    
    Row --> Split
    Split --> ExtractFields
    ExtractFields --> CreateID
    CreateID --> EntryList
    CreateID --> AFDict
```

**Sources:** [run_disobind.py:212-255]()

## UniProt Sequence Download

The system downloads UniProt sequences for all unique protein IDs in the input CSV. This process is parallelized using Python's `multiprocessing.Pool` for efficiency [run_disobind.py:261-278]().

### Sequence Download Process

1. **Extract Unique IDs:** The `get_unique_uni_ids()` method extracts all unique UniProt IDs from the entry list [run_disobind.py:280-295]().
2. **Parallel Download:** The `download_uniprot_seq()` method uses a multiprocessing pool to download sequences in parallel [run_disobind.py:299-327]().
3. **API Logic:** The `get_uniprot_seq()` function in `dataset/from_APIs_with_love.py` handles the HTTP requests to UniProt [dataset/from_APIs_with_love.py:10]().
4. **Caching:** Downloaded sequences are saved to `UniProt_seq.json` in the output directory to avoid re-downloading [run_disobind.py:108]().

**Sources:** [run_disobind.py:261-327](), [dataset/from_APIs_with_love.py:10-40]()

## FASTA File Generation

The `create_fasta_from_headers()` method generates FASTA files containing sequences for embedding generation [dataset/create_input_embeddings.py:186-235]().

### Global Scope (Default)

For global embeddings, the complete UniProt sequence is written to the FASTA file:

```
>P04273
MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGK...
```

**Sources:** [dataset/create_input_embeddings.py:207-210](), [run_disobind.py:66]()

## Embedding Generation

Disobind uses ProtT5 embeddings by default [run_disobind.py:64](). The embedding generation process involves calling external scripts from the ProtTrans library.

### Embedding Pipeline

```mermaid
flowchart TD
    FastaInput["FASTA File<br/>p1_p2_test.fasta"]
    
    GetEmb["get_embeddings()<br/>[utility.py:155]"]
    CheckType{"embedding_type?"}
    
    T5["ProtT5_embeddings()<br/>[utility.py:342]"]
    
    ChangeDir["os.chdir(ProtTrans/Embedding/)"]
    CallScript["subprocess.call(['python',<br/>'prott5_embedder.py',<br/>'--input', fasta_file,<br/>'--output', emb_file])"]
    
    EmbFile["Embeddings HDF5 File<br/>p1_p2_test.h5<br/>Per-residue embeddings<br/>Shape: [L, 1024]"]
    
    LoadH5["h5py.File() Open<br/>get_global_embeddings()<br/>[create_input_embeddings.py:276]"]
    SliceEmb["Slice Embeddings<br/>emb[start-1:end]<br/>for Fragment"]
    
    EmbDicts["prot1_emb Dictionary<br/>prot2_emb Dictionary<br/>{entry_id: np.array}"]
    
    FastaInput --> GetEmb
    GetEmb --> CheckType
    CheckType -->|T5| T5
    
    T5 --> ChangeDir
    ChangeDir --> CallScript
    CallScript --> EmbFile
    
    EmbFile --> LoadH5
    LoadH5 --> SliceEmb
    SliceEmb --> EmbDicts
```

### Embedding Properties

| Property | Value | Description |
|----------|-------|-------------|
| Embedding Type | ProtT5 | Default, set in `Disobind` class [run_disobind.py:64]() |
| Dimensions | [L, 1024] | Per-residue embeddings, L = sequence length |
| Storage Format | HDF5 | Key = entry ID or UniProt ID, Value = embedding array |
| Scope | Global | Uses full sequence for context [run_disobind.py:66]() |

**Sources:** [dataset/create_input_embeddings.py:239-273](), [dataset/utility.py:342-372](), [run_disobind.py:64-66]()

## Input Tensor Preparation

After embeddings are loaded, they must be converted to model-ready tensors. This involves padding, coarse-graining, and creating masks.

### Tensor Preparation Steps

The `get_input_tensors()` method [run_disobind.py:475-529]() performs the following operations:

1. **Calculate Effective Length:** Compute post-coarse-graining dimensions [run_disobind.py:488-490]().
2. **Create Target Mask:** Initialize a binary mask of ones with shape `[1, num_res1, num_res2]` [run_disobind.py:492]().
3. **Convert to Tensors:** Convert numpy arrays to PyTorch tensors and add batch dimension [run_disobind.py:504-507]().
4. **Apply Padding and Coarse-Graining:** Call `prepare_input()` from `src/utils.py` to pad sequences to 100 residues and apply coarse-graining [run_disobind.py:510-525]().
5. **Move to Device:** Transfer tensors to specified device (CPU or CUDA) [run_disobind.py:527]().

**Sources:** [run_disobind.py:475-529](), [src/utils.py:10-40]()

## Batch Processing

Predictions are processed in batches to manage memory efficiently. The default batch size is 200 protein pairs [run_disobind.py:78]().

**Diagram: Batch Processing Loop**

```mermaid
flowchart TD
    Start["Start: get_predictions()<br/>[run_disobind.py:168]"]
    CalcBatch["Calculate Batch Ranges<br/>np.arange(0, total_pairs, batch_size)"]
    
    Loop{"More<br/>Batches?"}
    SelectBatch["Select Batch:<br/>batch = prot_pairs[start:end]"]
    CreateEmb["create_embeddings(batch)<br/>[run_disobind.py:375]"]
    
    Predict["predict(required_tasks, af_dict)<br/>[run_disobind.py:532]"]
    ProcessEntry["For Each Entry in Batch:<br/>- Load embeddings<br/>- Get input tensors<br/>- Run model inference"]
    
    UpdateDict["Update predictions Dictionary"]
    Cleanup["Remove Processed Files:<br/>rm emb_file, fasta_file"]
    
    SaveOutput["Save predictions.npy"]
    End["End"]
    
    Start --> CalcBatch
    CalcBatch --> Loop
    Loop -->|Yes| SelectBatch
    SelectBatch --> CreateEmb
    CreateEmb --> Predict
    Predict --> ProcessEntry
    ProcessEntry --> UpdateDict
    UpdateDict --> Cleanup
    Cleanup --> Loop
    Loop -->|No| SaveOutput
    SaveOutput --> End
```

**Sources:** [run_disobind.py:168-207](), [run_disobind.py:532-661]()

## Error Handling and Validation

The input processing pipeline includes several validation steps:

1. **CSV Format Validation** [run_disobind.py:228-242](): Verifies each row has exactly 6 or 12 fields.
2. **UniProt Download Validation** [run_disobind.py:320-324](): Checks if sequence length is non-zero.
3. **Sequence Availability** [run_disobind.py:362-366](): Skips pairs with missing UniProt sequences and logs a warning.
4. **Coarse-Graining Validation** [run_disobind.py:144-151](): Verifies CG value is in [0, 1, 5, 10] and warns about potential terminal residue loss.

**Sources:** [run_disobind.py:144-151](), [run_disobind.py:228-242](), [run_disobind.py:320-324](), [run_disobind.py:362-366]()

---