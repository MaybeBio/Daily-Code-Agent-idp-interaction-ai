# OpenDNA desktop app

Tauri 2 shell with a Leptos (WebAssembly) frontend. The UI calls the Rust
library crates in `../crates` through typed Tauri commands. There is no
HTTP server and no open port.

```
app/
  protocol/   IPC message types shared by backend and UI
  desktop/    Tauri backend: commands, file import, report building
  ui/         Leptos UI, built with Trunk
  xtask/      developer tasks (fetching the Mol* viewer bundle)
```

## Prerequisites

- Rust 1.95+ with the `wasm32-unknown-unknown` target
- `cargo install trunk --locked` and `cargo install tauri-cli --version "^2" --locked`
- Linux only: WebKitGTK 4.1 and its build dependencies (see the Tauri prerequisites guide)

## Run

```bash
cargo run -p xtask -- fetch-molstar   # once: downloads Mol* 5.8.0, verifies its npm sha512
cd desktop && cargo tauri dev
```

`cargo tauri build` (from `desktop/`) produces platform installers.

## What works

- Import PDB, mmCIF and FASTA files. Each polymer chain becomes one entry in a local SQLite database under the OS app-data directory.
- Mol* 3D view. Hovering a residue in the sequence ruler highlights it; clicking focuses the camera on it.
- Sequence properties: ProtParam-compatible values and a Kyte–Doolittle profile.
- Structure analysis per chain: DSSP, relative solvent accessibility, φ/ψ, radius of gyration, total SASA, disulfide bonds and salt bridges.

Prediction and design tools arrive with `opendna-ml` and are not in the UI yet.
