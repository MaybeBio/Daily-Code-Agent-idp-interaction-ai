FROM alpine:3.18

LABEL Author="Luc Elliott" \
      Version="1.0" \
      Description="Lightweight HPC-friendly container with Linux utilities"

ENV PATH=/usr/local/bin:$PATH \
    LANG=C.UTF-8

# Install core utilities via apk (no static binary workarounds needed)
RUN apk add --no-cache \
    bash \
    busybox \
    wget \
    coreutils \
    grep \
    sed \
    tar \
    gzip

ENTRYPOINT ["/bin/sh", "-c"]