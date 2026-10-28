"""
Functional classification + alteration-geometry annotation for Pfam domain hits.

Replaces the previous kinase-only `is_kinase` flag with two orthogonal axes that
are more useful for hypothesis generation across the AS-isoform corpus:

  functional_class   what the domain DOES
                     {catalytic_kinase, catalytic_phosphatase,
                      catalytic_protease, catalytic_other, ub_ligase,
                      dna_binding, rna_binding, signaling_adapter,
                      oligomerization, other}

  alteration_type    WHERE the domain has changed (geometry of the alteration)
                     {intact, truncated_n, truncated_c, truncated_both,
                      internal_disruption, lost, gained}

Curation caveat
---------------
The Pfam clan/accession → functional_class mapping below is a curated
approximation, not a complete ontology.  It is **conservative**: a domain
gets a specific class only when the assignment is unambiguous.  Domains
landing in `other` are *not* less interesting — many transcription
co-factors, RNA-processing factors, and structural scaffolds end up here.

Two design notes:

  1. The same alteration_type can have *opposite* functional consequences
     depending on the domain — e.g. truncated_n of a SRC-family kinase
     removes the SH2/SH3 auto-inhibition (activating), but truncated_n
     of a typical receptor tyrosine kinase removes the catalytic lobe N
     (kinase-dead).  We therefore do NOT infer a direction (activating /
     dead / dominant-negative); we surface the geometry and let the
     researcher interpret.

  2. Multi-function domains (e.g. SAM is both oligomerization and
     adapter; SET is catalytic + binding) are assigned to their primary
     activity here.  Surfacing multiple classes per domain is left for a
     future iteration.
"""

from typing import Optional

# ── Curated clan → functional_class table ─────────────────────────────────────
# Each entry is annotated with the Pfam clan's common name for grep-ability.
# Only confidently-classified clans are listed; everything else falls through
# to ACCESSION_TO_CLASS, then to "other".
CLAN_TO_CLASS: dict[str, str] = {
    # Catalytic — kinases
    "CL0016": "catalytic_kinase",        # Pkinase superfamily (Ser/Thr + Tyr)

    # Catalytic — phosphatases
    "CL0031": "catalytic_phosphatase",   # Cysteine-based PTP superfamily
    "CL0064": "catalytic_phosphatase",   # His-phosphatase (HAD-like)

    # Catalytic — proteases
    "CL0035": "catalytic_protease",      # Papain-like cysteine proteases (CA clan)
    "CL0040": "catalytic_protease",      # Trypsin/chymotrypsin-like serine proteases
    "CL0123": "dna_binding",             # Homeodomain   (assigned below, kept here as comment)

    # Catalytic — other enzymatic
    "CL0010": "catalytic_other",         # Rossmann fold dehydrogenases
    "CL0023": "catalytic_other",         # P-loop NTPase (helicases, ATPases, GTPases)
    "CL0063": "catalytic_other",         # FAD/NAD binding domain (often catalytic)

    # DNA-binding
    "CL0072": "dna_binding",             # Helix-turn-helix superfamily
    "CL0114": "dna_binding",             # Zinc finger superclan
    "CL0481": "dna_binding",             # bZIP-like

    # RNA-binding
    "CL0221": "rna_binding",             # RRM-like
    "CL0099": "rna_binding",             # KH-like
    "CL0036": "rna_binding",             # Sm-like / LSm (snRNP component)

    # Signaling adapters / scaffolds
    "CL0254": "signaling_adapter",       # SH2-like
    "CL0010" : "catalytic_other",        # (already set) Rossmann fold dehydrogenases

    # Oligomerization / repeats / scaffold
    "CL0020": "oligomerization",         # TPR-like repeats
    "CL0465": "oligomerization",         # Ankyrin repeat
}

# Homeodomain has its own clan in some Pfam releases; assign explicitly here so
# CL0123 lands in dna_binding regardless of the upstream clan curation drift.
CLAN_TO_CLASS["CL0123"] = "dna_binding"


# ── Accession-level overrides / additions ──────────────────────────────────────
# Used for families whose clan is missing, ambiguous, or absent in our DB's
# pfam_clan column (Pfam-A.hmm.dat parsing fails for some families).
ACCESSION_TO_CLASS: dict[str, str] = {
    # Kinases (defensive — these are covered by CL0016 but listed for safety)
    "PF00069": "catalytic_kinase",   # Pkinase (Ser/Thr)
    "PF07714": "catalytic_kinase",   # Pkinase_Tyr
    "PF01163": "catalytic_kinase",   # RIO1 atypical kinase
    "PF02816": "catalytic_kinase",   # Alpha-kinase

    # Phosphatases
    "PF00102": "catalytic_phosphatase",   # Y_phosphatase (PTP)
    "PF00782": "catalytic_phosphatase",   # DSPc (dual-specificity)
    "PF00481": "catalytic_phosphatase",   # PP2C
    "PF00149": "catalytic_phosphatase",   # Metallophos (PP1/PP2A catalytic)

    # Proteases
    "PF00089": "catalytic_protease",   # Trypsin
    "PF00112": "catalytic_protease",   # Peptidase_C1 (papain)
    "PF00656": "catalytic_protease",   # Caspase
    "PF00026": "catalytic_protease",   # Asp_protease (aspartyl/pepsin)
    "PF00207": "catalytic_protease",   # Alpha-2-macroglobulin (broad protease inhibitor — keep tagged)
    "PF00227": "catalytic_protease",   # Proteasome A/B subunits

    # Ubiquitin/SUMO machinery
    "PF00097": "ub_ligase",   # zf-C3HC4 (RING)
    "PF13639": "ub_ligase",   # zf-RING_2
    "PF00179": "ub_ligase",   # UQ_con (E2)
    "PF00444": "ub_ligase",   # HECTc (HECT catalytic)
    "PF00632": "ub_ligase",   # HECT
    "PF00917": "ub_ligase",   # MATH/TRAF
    "PF00240": "other",       # Ubiquitin (the moiety itself — substrate, not enzyme)

    # DNA-binding
    "PF00046": "dna_binding",   # Homeobox
    "PF00010": "dna_binding",   # HLH
    "PF00170": "dna_binding",   # bZIP_1
    "PF07716": "dna_binding",   # bZIP_2
    "PF00096": "dna_binding",   # C2H2 zinc finger
    "PF13912": "dna_binding",   # C2H2-type zinc finger
    "PF00104": "dna_binding",   # Hormone receptor LBD (also ligand-binding, but DBD context)
    "PF00105": "dna_binding",   # zf-C4 (nuclear receptor)
    "PF00250": "dna_binding",   # Forkhead
    "PF00292": "dna_binding",   # PAX (paired box)
    "PF00352": "dna_binding",   # TFIID/TBP

    # RNA-binding
    "PF00076": "rna_binding",   # RRM_1
    "PF13893": "rna_binding",   # RRM_5
    "PF00013": "rna_binding",   # KH_1
    "PF07650": "rna_binding",   # KH_2
    "PF00035": "rna_binding",   # dsrm (dsRNA-binding)
    "PF00098": "rna_binding",   # zf-CCHC (RNA-binding zinc knuckle)

    # RNA helicases / processing (catalytic but routed via RNA pathway)
    "PF00270": "catalytic_other",   # DEAD helicase
    "PF00271": "catalytic_other",   # Helicase_C

    # Signaling adapters / scaffolds (not enzymatic, mediate interactions)
    "PF00017": "signaling_adapter",   # SH2
    "PF00018": "signaling_adapter",   # SH3_1
    "PF14604": "signaling_adapter",   # SH3_9
    "PF00595": "signaling_adapter",   # PDZ
    "PF00169": "signaling_adapter",   # PH
    "PF00640": "signaling_adapter",   # PID / PTB
    "PF00397": "signaling_adapter",   # WW
    "PF00498": "signaling_adapter",   # FHA
    "PF00536": "signaling_adapter",   # SAM_1 (also dimerization)

    # Oligomerization / structural scaffolds
    "PF00023": "oligomerization",   # Ankyrin (also CL0465)
    "PF12796": "oligomerization",   # Ank_2
    "PF13637": "oligomerization",   # Ank_4
    "PF00400": "oligomerization",   # WD40
    "PF00643": "oligomerization",   # B-box zinc finger (TRIM family scaffold)
    "PF02184": "oligomerization",   # HAT (TPR-related)
    "PF13176": "oligomerization",   # TPR_7
    "PF12895": "oligomerization",   # Ankyrin_rpt-containing
    "PF01400": "oligomerization",   # Astacin (also protease — kept here for now)
}

# Subset that *was* the old is_kinase=True set, to preserve legacy semantics.
LEGACY_KINASE_ACCESSIONS = {"PF00069", "PF07714", "PF00081"}
LEGACY_KINASE_CLANS = {"CL0016"}


# ── Public API ────────────────────────────────────────────────────────────────

def classify_domain(pfam_accession: Optional[str],
                    pfam_clan: Optional[str]) -> str:
    """
    Return the functional_class label for a Pfam hit.

    Resolution order:
      1. accession-level override (most specific)
      2. clan-level mapping
      3. `other` (default)

    Returns one of the values listed in this module's docstring.
    """
    if pfam_accession:
        # Strip any version suffix; our DB stores without version, but be defensive
        acc = pfam_accession.split(".")[0]
        cls = ACCESSION_TO_CLASS.get(acc)
        if cls is not None:
            return cls
    if pfam_clan:
        cls = CLAN_TO_CLASS.get(pfam_clan)
        if cls is not None:
            return cls
    return "other"


def is_catalytic_class(functional_class: str) -> bool:
    """True if the class is enzymatic (covers catalytic_* + ub_ligase)."""
    return functional_class.startswith("catalytic_") or functional_class == "ub_ligase"


# ── Alteration geometry ───────────────────────────────────────────────────────
#
# These thresholds describe what we count as N- vs C- vs internal disruption
# based on the HMM model coordinates of the highest-scoring hit:
#
#     0%                                                              100%
#     ├─────────────────────── HMM model ──────────────────────────────┤
#                  ↑                                  ↑
#               hmm_start                          hmm_end
#
# We say the N-terminal portion is "missing" if hmm_start is well past the
# beginning of the HMM (i.e. our query sequence aligned to a partial model
# starting well into it), and similarly for the C-terminal end.

N_TERM_PAD = 0.10   # if hmm_start > 10% of HMM, call N-terminal portion missing
C_TERM_PAD = 0.90   # if hmm_end   < 90% of HMM, call C-terminal portion missing


def alteration_type(domain_status: str,
                    hmm_start: Optional[int],
                    hmm_end: Optional[int],
                    hmm_length: Optional[int],
                    domain_coverage: Optional[float]) -> str:
    """
    Derive the alteration-geometry label from the per-row HMM coordinates +
    the existing `domain_status` classification.

    The mapping is:
      domain_status == 'present'    → 'intact'        (full or near-full coverage)
      domain_status == 'truncated'  → one of truncated_n / truncated_c /
                                       truncated_both / internal_disruption
      domain_status == 'lost'       → 'lost'
      domain_status == 'gained'     → 'gained'
      anything else                 → 'unknown'
    """
    if domain_status == "lost":
        return "lost"
    if domain_status == "gained":
        return "gained"
    if domain_status == "present":
        return "intact"
    if domain_status != "truncated":
        return "unknown"

    # Truncated: figure out *which* end(s) are missing
    if not hmm_length or hmm_length <= 0:
        return "truncated_unknown"

    n_cut = hmm_start is not None and hmm_start > hmm_length * N_TERM_PAD
    c_cut = hmm_end   is not None and hmm_end   < hmm_length * C_TERM_PAD

    if n_cut and c_cut:
        # Both ends inside the HMM model — but is this because the middle is
        # actually intact (real bilateral truncation) or because the hit covers
        # a thin internal slice?  Use coverage to disambiguate.
        if domain_coverage is not None and domain_coverage < 0.50:
            return "internal_disruption"
        return "truncated_both"
    if n_cut:
        return "truncated_n"
    if c_cut:
        return "truncated_c"
    # Neither end shifted — but coverage flagged truncation.  This generally
    # means an alignment gap inside the model: surface as internal.
    return "internal_disruption"
