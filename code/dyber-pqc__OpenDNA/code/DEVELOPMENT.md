# OpenDNA Development Guide

## Workspaces

Four Cargo workspaces, deliberately separate so that the application does not
depend on what the data and training tools need — Arrow, Parquet and an HTTP
client have no business in a desktop binary.

| Path | Workspace | Contents |
|---|---|---|
| `crates/` | root | core types, PDB/mmCIF/FASTA I/O, structure analysis, storage, ESM-2 and the prediction heads |
| `app/` | `app` | Tauri 2 desktop shell, Leptos UI, the IPC protocol between them |
| `datasets/` | `datasets` | `opendna-data`, which builds the datasets from pinned upstream sources |
| `training/` | `training` | embedding cache, head training, model cards, external benchmarks |
| `models/` | — | trained heads and their cards; weights are never committed |

## Build and test

```bash
cargo test --workspace
```

Each other workspace is built from its own directory:

```bash
cd datasets && cargo test --workspace
cd training && cargo test --workspace
```

The desktop application needs a WebAssembly target and Trunk:

```bash
cd app
rustup target add wasm32-unknown-unknown
cargo install trunk
cargo run -p xtask -- fetch-molstar
cd ui && trunk build --release && cd ..
cargo run -p opendna-desktop
```

GPU builds need the `cuda` feature and a toolkit version candle supports. On
Windows, `scripts/cargo-cuda.ps1` puts nvcc and MSVC on `PATH` and explains the
version bounds; pass `-j 4` on the first build, because one nvcc per core is
enough to exhaust memory on a many-core machine.

## Conventions

- Minimum Rust 1.95, declared identically in every workspace and checked in CI.
  A declared minimum that nothing builds against is a guess.
- `unsafe_code = "forbid"` across the workspace. Among other things this rules
  out memory-mapped safetensors loading, which is why weights are read rather
  than mapped.
- `cargo fmt`, and `clippy` at `pedantic` with warnings denied. Where a lint is
  allowed, the `allow` carries a comment saying why.
- Errors are `Result`; a missing model or tool is `Err(Unavailable)`.

## Ground rules

These are the rules the rewrite exists to enforce, and they are not style
preferences:

- **No fake science.** If a model or tool is absent, the answer is an error
  saying so. Never a heuristic presented as a prediction. The Python layer this
  replaced fell back to hydropathy arithmetic and reported the result in
  kcal/mol, which is the specific failure being designed out.
- **Report the baseline.** Every model card states what the model is compared
  against. A model that does not beat its baseline is not published.
- **Test once.** Validation chooses the epoch; the test split is scored
  afterwards, once. Nothing is tuned against a test number, and nothing is
  tuned against an external benchmark at all.
- **Say what the data does not cover.** Every card's limitations section is
  written by hand and re-read before release.
- **Offline by default.** No open port, no service. The first prediction
  downloads a checkpoint; after that the application works with no network.
