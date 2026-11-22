# The OpenDNA command-line tool in a container.
#
# For a scheduler. Nextflow and Snakemake both take a container per process,
# and a pipeline that runs parts in parallel wants every part running the same
# build of the same binary — which is the whole reason `opendna merge` reads
# each part's manifest and refuses a mixture.
#
# It builds the CLI only. The desktop application is a GUI and does not belong
# in a container, and the previous version of this file built a Python API
# server from directories that no longer exist.
#
#   docker build -t opendna .
#   docker run --rm -v "$PWD:/work" -v /mnt/models:/models:ro opendna \
#       --models /models predict --task disorder -i /work/in.fasta -o /work/out.csv
#
# Nothing is baked in but the binary. No weights, no network service, no port.
#
# Built and smoke-tested by the `container` job in CI on every commit, because
# it could not be built on the machine it was written on: Docker is installed
# there but the daemon is not running. The CI job is the verification.

FROM rust:1.97-slim-bookworm AS builder

# A C toolchain and cmake for the dependencies that build native code. Not for
# OpenDNA itself, which forbids unsafe code and has no build script.
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential pkg-config cmake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY Cargo.toml Cargo.lock clippy.toml rustfmt.toml ./
COPY crates ./crates

# --locked, so the image is built from the dependency versions the repository
# pins rather than whatever resolves on the day.
RUN cargo build --release --locked -p opendna-cli

FROM debian:bookworm-slim

# Certificate roots, for the one thing the tool reaches the network for: the
# first download of the embedding model. Stage a bundle and pass --offline and
# it never needs them.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/target/release/opendna /usr/local/bin/opendna
COPY LICENSE /usr/share/doc/opendna/LICENSE
COPY crates/opendna-cli/README.md /usr/share/doc/opendna/README.md

# Not root: nothing here needs it, and a pipeline mounting a shared filesystem
# should not have a container writing to it as root.
RUN useradd --create-home --uid 1000 opendna
USER opendna
WORKDIR /work

# Where `--models` is expected to be mounted, and where the Hugging Face cache
# would otherwise land in a home directory that a scheduler may not persist.
ENV OPENDNA_MODELS=/models \
    HF_HOME=/work/.cache/huggingface

ENTRYPOINT ["opendna"]
CMD ["--help"]
