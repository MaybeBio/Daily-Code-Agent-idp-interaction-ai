# Rust crates

OpenDNA is being rewritten in Rust. These crates are the new core. The
Python backend in `python/` and the React UI in `ui/src/` stay in place until
the rewrite reaches feature parity, then they are removed.

| Crate | Purpose | Status |
|---|---|---|
| `opendna-core` | Amino acids, sequences, structures, geometry | usable |
| `opendna-io` | PDB read/write, mmCIF read, FASTA read/write | usable |
| `opendna-analysis` | ProtParam properties, φ/ψ/ω, DSSP, Shrake–Rupley SASA, Kabsch superposition, disulfides, salt bridges | usable |
| `opendna-store` | SQLite protein store, sequence version history | usable |
| `opendna-hal` | Hardware detection for model sizing | usable; CUDA and Apple only |

The desktop app built on these crates lives in [`app/`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/app/README.md).

Planned: `opendna-ml` (candle models), `opendna-engines`, `opendna-jobs`,
`opendna-provenance`, `opendna-workflow`, `opendna-auth` and `opendna-server`.

## Conventions

- No result without a method. If a quantity can't be computed (missing
  model, non-standard residue), return an error instead of an estimate.
- Analyses are tested against published values on real PDB entries
  (`crates/opendna-analysis/tests/`), not only synthetic inputs.
- `cargo fmt`, `cargo clippy --all-targets -- -D warnings` (pedantic),
  `cargo test` and `cargo deny check` must pass. CI runs all four on
  Linux, macOS and Windows.

## Known limitations

- DSSP: β-bulges aren't merged into ladders and PP-II helices aren't
  assigned. Everything else follows Kabsch & Sander (1983).
- The mmCIF reader reads `_atom_site` only. There is no mmCIF writer yet.
- Instability index is not implemented. The earlier Python version used
  an incomplete weight table, so it was dropped rather than ported.
