# Contributing to OpenDNA

Thanks for your interest! OpenDNA is community-driven. Here's how to help.

## Ways to contribute

- **Report bugs** — open an issue with reproduction steps
- **Request features** — open an issue describing the use case
- **Improve documentation** — typos, clarifications, new tutorials
- **Add tests** — coverage on engines is the highest-impact contribution
- **Build features** — pick something from the roadmap or propose your own
- **Share workflows** — submit example pipelines to the cookbook
- **Validate science** — run benchmarks and report results

## Quick start for developers

```bash
git clone https://github.com/dyber-pqc/OpenDNA.git
cd OpenDNA

# The libraries
cargo test --workspace

# The desktop application
cd app
rustup target add wasm32-unknown-unknown
cargo install trunk
cargo run -p xtask -- fetch-molstar
cd ui && trunk build --release && cd ..
cargo run -p opendna-desktop
```

There is no server to run and no Python or Node toolchain to install. See
[DEVELOPMENT.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/DEVELOPMENT.md) for the workspace layout, the GPU build, and
the ground rules a change is held to.

## Code conventions

- **Rust**: `cargo fmt`; `clippy` at `pedantic` with warnings denied. An
  `allow` carries a comment saying why.
- **Errors**: `Result` everywhere. A missing model or tool is an error that
  names what is missing — never a heuristic standing in for a prediction.
- **CSS**: use the tokens in `app/ui/styles/tokens.css`, no inline styles.
- **Commit messages**: say what changed and why it was wrong before. Not
  Conventional Commits prefixes.

## Pull request process

1. Fork the repo and create a branch from `main`
2. Make your changes with clear commits
3. Add tests for new functionality
4. Run the test suite and make sure it passes
5. Update documentation if user-facing
6. Open a PR with a clear description
7. Wait for CI to pass
8. Address review feedback
9. Maintainer merges

## Things we will NOT accept

- Code without tests for new functionality
- Breaking changes without discussion
- Proprietary code we can't relicense
- Vendoring of large binaries
- Changes that break cross-platform support
- Anything with hardcoded secrets
- AI-generated code without verification (we love AI tools, but you must read and verify the output)

## Code of Conduct

By participating, you agree to follow our [Code of Conduct](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/CODE_OF_CONDUCT.md).

## Licence

By contributing, you agree to licence your contributions under the same terms
as the project: **[Apache-2.0](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/LICENSE)**.

This previously said "Apache 2.0 + Commons Clause", which is a different and
more restrictive grant and is not the licence in the `LICENSE` file or in any
crate's metadata. It mattered more here than in a README: this paragraph is
what a contributor agrees to, and it was asking them to agree to terms the
project does not use.

Model weights are a separate question — they carry the licence their training
data requires, and those differ. See the [licence section of the
README](README.md#licence).

## Questions?

- Open a GitHub Discussion
- File an issue with the `question` label
- Join our Discord (link coming soon)
