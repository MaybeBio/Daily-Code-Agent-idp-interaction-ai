#!/usr/bin/env bash
#
# Sign a release's SHA256SUMS with GPG, producing a detached signature.
#
# This is the cheapest useful signature and the only one that covers every
# artefact on every platform at once: the checksum file names all of them, so
# one signature over it vouches for all of them. It costs nothing but a key.
#
# What it does not do is stop the operating system warning about an unsigned
# installer. That needs platform code signing -- Authenticode on Windows, a
# Developer ID and notarization on macOS -- which needs certificates that cost
# money. docs/SIGNING.md says which, what they cost, and what each one buys.
#
#   scripts/sign-release.sh release/        # sign release/SHA256SUMS
#   scripts/sign-release.sh --verify release/
#
# Signing reads the private key from GPG_SIGNING_KEY, an ASCII-armoured secret
# key, and the passphrase from GPG_PASSPHRASE if the key has one. The key is
# imported into a temporary keyring that is deleted afterwards, so this never
# touches the keyring of whoever runs it and never writes the key to a file
# that outlives the run.
#
# Verifying needs no secret: it uses the public key committed at
# docs/signing-key.asc, so anybody can run the verify path.

set -euo pipefail

mode=sign
if [ "${1:-}" = "--verify" ]; then
    mode=verify
    shift
fi

dir="${1:-release}"
sums="$dir/SHA256SUMS"
signature="$sums.asc"
# Overridable so the round-trip test in CI can verify against a key it
# generated, rather than needing the real one.
public_key="${OPENDNA_SIGNING_PUBKEY:-$(dirname "$0")/../docs/signing-key.asc}"

if [ ! -f "$sums" ]; then
    echo "no $sums. Run this against the directory holding the release artefacts." >&2
    exit 2
fi

# A keyring of its own, removed on the way out however this exits.
home="$(mktemp -d)"
trap 'rm -rf "$home"' EXIT
chmod 700 "$home"
export GNUPGHOME="$home"

if [ "$mode" = verify ]; then
    if [ ! -f "$signature" ]; then
        echo "no $signature to check" >&2
        exit 2
    fi
    if [ ! -f "$public_key" ]; then
        echo "no public key at $public_key, so there is nothing to verify against." >&2
        echo "Commit the public half of the signing key there." >&2
        exit 2
    fi
    gpg --quiet --import "$public_key"
    # --status-fd so the result is read from gpg's machine-readable output
    # rather than from prose that changes between versions.
    if gpg --verify --status-fd 1 "$signature" "$sums" 2>/dev/null | grep -q '^\[GNUPG:\] VALIDSIG'; then
        echo "SHA256SUMS signature is good."
        # And the checksums themselves, since a good signature over wrong
        # checksums is not what anyone wanted to know.
        (cd "$dir" && sha256sum --check --ignore-missing SHA256SUMS)
    else
        echo "SHA256SUMS signature did NOT verify. Do not use these files." >&2
        exit 1
    fi
    exit 0
fi

if [ -z "${GPG_SIGNING_KEY:-}" ]; then
    echo "GPG_SIGNING_KEY is not set, so nothing was signed." >&2
    echo "This is not a failure on a build that is not meant to be signed." >&2
    exit 3
fi

printf '%s' "$GPG_SIGNING_KEY" | gpg --quiet --batch --import
fingerprint="$(gpg --list-secret-keys --with-colons | awk -F: '/^fpr:/ { print $10; exit }')"
if [ -z "$fingerprint" ]; then
    echo "GPG_SIGNING_KEY imported no secret key." >&2
    exit 1
fi

rm -f "$signature"
gpg --batch --yes --pinentry-mode loopback \
    ${GPG_PASSPHRASE:+--passphrase "$GPG_PASSPHRASE"} \
    --local-user "$fingerprint" \
    --armor --detach-sign --output "$signature" "$sums"

# Check our own output before claiming success. A signature nobody verified is
# a file, not a signature.
if gpg --verify --status-fd 1 "$signature" "$sums" 2>/dev/null | grep -q '^\[GNUPG:\] VALIDSIG'; then
    echo "Signed $sums with $fingerprint"
    echo "Wrote $signature"
else
    echo "the signature this just wrote does not verify" >&2
    exit 1
fi
