#!/bin/bash

############################## INPUTS ##############################
# $1 --> input PDB file
# $2 --> ligand resname                (optional)
# $3 --> ligand charge                 (optional)
# $4 --> ph for ligand protonation (normally 7.4)
#
# Options (environment variables, see README):
#   MODEL=ff19SB_tip3p     ff19SB_tip3p | a99SBdisp_tip4pd
#   WATER=authors          a99SBdisp_tip4pd only: authors | amber
#   BOX=rect               rect | cubic
#   PAD=15.0               solvent padding (A)
#   MINBOX=<none>          target MINIMUM box dimension (A). Overrides PAD with
#                          per-axis padding so the box is >= MINBOX on every
#                          axis - for an elongated solute that tumbles, the
#                          limit is the SHORTEST box side vs the solute's
#                          DIAMETER, not the padding. Needs BOX=rect.
#   PADFLOOR=8.0           with MINBOX: smallest per-axis padding allowed
#   SS_CUTOFF=2.5          SG-SG distance (A) below which CYS pairs are bonded
#   CN_CUTOFF=2.0          C(i)-N(i+1) distance (A) above which the backbone is
#                          treated as broken (TER inserted, no peptide bond)
#   SALT=<model default>   salt (M) on top of neutralization: ff19SB_tip3p 0, a99SBdisp_tip4pd 0.15
#   CATION=<model default> Na+ (ff19SB_tip3p) | K+ (a99SBdisp_tip4pd)
#
# Example: MODEL=a99SBdisp_tip4pd BOX=cubic bash utils/build.sh input.pdb

PDBFILE="$1"
RESNAME="$2"
LIGCHARGE="$3"
LIGPH="$4"

MODEL=${MODEL:-ff19SB_tip3p}
BOX=${BOX:-rect}
PAD=${PAD:-15.0}

set -e

# Check if at least the PDB file is provided
if [ -z "$PDBFILE" ]; then
    echo "Usage: [MODEL=...] $0 <pdb_file> [ligand_resname ligand_charge ph]"
    exit 1
fi

# Model settings
case "$MODEL" in
    ff19SB_tip3p)
        WATER=tip3p
        LEAPRC=leaprc.pipe_ff19SB_tip3p;   WATERBOX=TIP3PBOX
        SALT=${SALT:-0};    CATION=${CATION:-Na+};  ADD8=false ;;
    ff19SB_tip4pew)
        # Native TIP4P-Ew (NOT the TIP4P-D of a99SBdisp_tip4pd's WATER=authors, which only
        # reuses the TIP4PEWBOX starting box name but loads different, D-variant parameters).
        # ff19SB was parameterised/validated against OPC, not TIP4P-Ew; added here anyway
        # because other simulations in this project standardise on TIP4P-Ew.
        WATER=tip4pew
        LEAPRC=leaprc.pipe_ff19SB_tip4pew; WATERBOX=TIP4PEWBOX
        SALT=${SALT:-0.15}; CATION=${CATION:-Na+};  ADD8=false ;;
    a99SBdisp_tip4pd)
        WATER=${WATER:-authors}
        case "$WATER" in
            authors) WATERBOX=TIP4PEWBOX ;;
            amber)   WATERBOX=TIP4PDBOX ;;
            *) echo "ERROR: WATER must be 'authors' or 'amber'"; exit 1 ;;
        esac
        LEAPRC=leaprc.pipe_a99SBdisp_tip4pd_${WATER}
        SALT=${SALT:-0.15}; CATION=${CATION:-K+};   ADD8=true ;;
    *) echo "ERROR: MODEL must be 'ff19SB_tip3p', 'ff19SB_tip4pew' or 'a99SBdisp_tip4pd'"; exit 1 ;;
esac
if [ "$BOX" = "cubic" ]; then ISO=iso; elif [ "$BOX" = "rect" ]; then ISO=""; else echo "ERROR: BOX must be 'rect' or 'cubic'"; exit 1; fi

# Create main directory to save the preparation steps
FOLDER=build
mkdir -p $FOLDER

echo "MODEL=$MODEL WATER=$WATER BOX=$BOX PAD=$PAD SALT=$SALT CATION=$CATION" | tee $FOLDER/build_settings.txt

# Our force-field files all start with pipe_ : tleap searches its own folders BEFORE -I,
# so a file with an AmberTools name would be silently replaced by the AmberTools one
LEAP="tleap -s -I utils/forcefields/$MODEL"
fill () {  # fill the placeholders of a leap template
    sed -e "s/__LEAPRC__/$LEAPRC/" -e "s/__WATERBOX__/$WATERBOX/" -e "s/__PAD__/$PAD/" -e "s/__ISO__/$ISO/" \
        -e "s/__CATION__/$CATION/g" -e "s/__SALT__/$SALT/" -e "s/__NSALT__/$NSALT/g" -e "s/__SALTLINE__/$SALTLINE/" "$1" \
    | sed -e "/__DISULFIDES__/r $FOLDER/disulfides.leap" -e "/__DISULFIDES__/d"
}

# Flag to indicate if ligand is present
LIGAND_PRESENT=false

# If all three inputs are provided, process the ligand
if [ -n "$RESNAME" ] && [ -n "$LIGCHARGE" ]; then
    LIGAND_PRESENT=true
    echo "Ligand detected: $RESNAME with charge $LIGCHARGE"

    # Extract the ligand only from the PDB
    grep -w "$RESNAME" "$PDBFILE" | sed "s/$RESNAME/out/g" > "$FOLDER/out_noh.pdb"

    # Protonate the ligand
    python utils/protonated_pdb.py "$FOLDER"/out_noh.pdb "$FOLDER"/out_ph.pdb $LIGPH
    sed -i 's/UNL/out/g' "$FOLDER"/out_ph.pdb

    # Run antechamber to generate ligand parameters
    antechamber -i "$FOLDER"/out_ph.pdb -fi pdb -o "$FOLDER"/out.mol2 -fo mol2 -c bcc -rn out -at gaff2 -nc "$LIGCHARGE"
    rm -f ANTECHAMBER_* ATOMTYPE.INF
    mv sqm* "$FOLDER"

    # Generate PREPI file
    antechamber -i "$FOLDER"/out_ph.pdb -fi pdb -o "$FOLDER"/out.prepin -fo prepi -c bcc -at gaff2 -nc "$LIGCHARGE"
    rm -f ANTECHAMBER_* ATOMTYPE.INF
    mv sqm* "$FOLDER"

    # Generate FRCMOD file
    parmchk2 -i "$FOLDER"/out.mol2 -f mol2 -o "$FOLDER"/out.frcmod -s gaff2
    rm -f PREP.INF NEWPDB.PDB

    # Convert the mol2 file back to PDB
    # to ensure compatibility between antechamber atom names and the PDB file
    # before running tleap
    # -dr no: pure format conversion. antechamber writes every bond order as 1
    # in out.mol2, so acdoctor would reject any C=O / carboxylate as an
    # "open valence" here, although the atom types and charges are correct.
    antechamber -i "$FOLDER"/out.mol2 -fi mol2 -o "$FOLDER"/out.pdb -fo pdb -dr no

    LEAPFILE=utils/build_complex.leap
else
    echo "No ligand detected. Only protein will be processed."
    LEAPFILE=utils/build_protein.leap
fi

# Extract protein without hydrogens
vmd -dispdev text -e utils/extract_protein.tcl -args "$PDBFILE" "$FOLDER/protein_noh.pdb" > $FOLDER/extract_protein.log

# VMD's writepdb drops all TER records, and tleap bonds across every break it
# is not told about (chain ends AND backbone gaps). Put TERs back at chain-ID
# changes and at real gaps (C(i)-N(i+1) > 2.0 A, overridable with CN_CUTOFF=...).
python utils/mark_chain_breaks.py "$FOLDER/protein_noh.pdb" "${CN_CUTOFF:-2.0}" | tee $FOLDER/chain_breaks.log

# Detect and bond disulfides (SG-SG < 2.5 A, overridable with SS_CUTOFF=...).
# Renames the matching CYS -> CYX in $FOLDER/protein_noh.pdb (in place) and
# writes $FOLDER/disulfides.{leap,json}; both are read above by fill().
python utils/find_disulfides.py "$FOLDER/protein_noh.pdb" "$FOLDER" "${SS_CUTOFF:-2.5}" | tee $FOLDER/disulfides.log

# ---- MINBOX: per-axis padding to a target minimum box dimension ----------
# solvatebox's box side is linear in the padding: side_i = vdw_extent_i + 2*pad.
# One cheap solvate pass at a known padding therefore measures vdw_extent_i
# exactly, and the padding that puts every side at MINBOX follows. This beats
# a hard-coded slack constant, which would drift with solvent/closeness.
if [ -n "${MINBOX:-}" ]; then
    if [ "$BOX" != "rect" ]; then
        echo "ERROR: MINBOX needs BOX=rect (it sets the per-axis padding itself)"; exit 1
    fi
    CALPAD=5.0
    echo "Calibrating per-axis padding for MINBOX=$MINBOX A (probe pad $CALPAD)"
    PAD_SAVED="$PAD"; PAD="$CALPAD"
    fill $LEAPFILE | sed '/#---- IONS ----/,$d' > $FOLDER/build_calib.leap
    echo "quit" >> $FOLDER/build_calib.leap
    $LEAP -f $FOLDER/build_calib.leap > $FOLDER/tleap_calib.log
    rm -f leap.log
    CALBOX=$(grep -oP 'Total vdw box size:\s+\K[0-9.]+\s+[0-9.]+\s+[0-9.]+' \
             $FOLDER/tleap_calib.log | tail -1)
    if [ -z "$CALBOX" ]; then
        echo "ERROR: no box size in $FOLDER/tleap_calib.log"; exit 1
    fi
    PAD=$(python - "$CALBOX" "$CALPAD" "$MINBOX" "${PADFLOOR:-8.0}" <<'PY'
import sys
box = [float(v) for v in sys.argv[1].split()]
cal, want, floor = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
ext = [b - 2 * cal for b in box]                 # solute vdw extent per axis
pads = [max(floor, (want - e) / 2.0) for e in ext]
sides = [e + 2 * p for e, p in zip(ext, pads)]
print("{%s}" % " ".join(f"{p:.2f}" for p in pads))
print("  solute vdw extent : %s A" % "  ".join(f"{e:.1f}" for e in ext),
      file=sys.stderr)
print("  per-axis padding  : %s A" % "  ".join(f"{p:.1f}" for p in pads),
      file=sys.stderr)
print("  resulting box     : %s A  (target >= %.1f)"
      % ("  ".join(f"{s:.1f}" for s in sides), want), file=sys.stderr)
PY
)
    echo "MINBOX=$MINBOX PADFLOOR=${PADFLOOR:-8.0} PAD=$PAD (was PAD=$PAD_SAVED)" \
        | tee -a $FOLDER/build_settings.txt
fi

# System preparation using tleap
NSALT=0; SALTLINE="# "
if python -c "import sys; sys.exit(0 if float('$SALT') > 0 else 1)"; then
    # Pass 1: solvate only, to count the waters for the salt concentration
    echo "Running tleap pass 1 (water count for salt)"
    fill $LEAPFILE | sed '/#---- IONS ----/,$d' > $FOLDER/build_pass1.leap
    echo "quit" >> $FOLDER/build_pass1.leap
    $LEAP -f $FOLDER/build_pass1.leap > $FOLDER/tleap_pass1.log
    NWAT=$(grep -oP 'Added \K[0-9]+(?= residues)' $FOLDER/tleap_pass1.log | tail -1)
    if [ -z "$NWAT" ]; then echo "ERROR: could not read water count from $FOLDER/tleap_pass1.log"; exit 1; fi
    # ion pairs = SALT [M] * N_water / 55.51 [M]
    NSALT=$(python -c "print(int(round($SALT * $NWAT / 55.51)))"); SALTLINE=""
    echo "Waters: $NWAT  ->  ${CATION}/Cl- pairs for ${SALT} M: $NSALT"
fi

echo "Running tleap ($LEAPFILE)"
fill $LEAPFILE > $FOLDER/build_system.leap
$LEAP -f $FOLDER/build_system.leap > $FOLDER/tleap_build.log
mv leap.log "$FOLDER" 2>/dev/null || true
if grep -q "Errors = [1-9]" $FOLDER/tleap_build.log || [ ! -s step3_input.parm7 ]; then
    echo "ERROR: tleap failed, see $FOLDER/tleap_build.log"; exit 1
fi
# Every file tleap loaded, with its full path (provenance)
grep -E "^(Loading parameters|Loading library|----- Source):" $FOLDER/tleap_build.log | sort -u > $FOLDER/loaded_files.txt

# a99SB-disp: add the extra backbone dihedrals (as in pipeline_amber_tip4pd's add8.py)
if [ "$ADD8" = true ]; then
    mv step3_input.parm7 $FOLDER/step3_input_noadd8.parm7
    python utils/add8.py $FOLDER/step3_input_noadd8.parm7 step3_input.parm7 > $FOLDER/add8.log
fi

# Final step: create restraints and move box to 0.0.0
mkdir -p restraints
vmd -dispdev text -e utils/create_omm.tcl > $FOLDER/create_omm.log

# Check that the intended parameters ended up in the topology
python utils/check_topology.py step3_input.parm7 $MODEL $WATER $FOLDER/disulfides.json | tee $FOLDER/topology_report.txt
if grep -q "^FAIL" $FOLDER/topology_report.txt; then
    echo "ERROR: topology check failed (see $FOLDER/topology_report.txt)"; exit 1
fi

# OpenMM inputs for this model (edit these copies, not the templates in inputs/)
rm -f step4.*_equilibration.inp step5_production.inp
cp inputs/$MODEL/*.inp .
echo "Copied OpenMM inputs from inputs/$MODEL/: $(ls step4.*_equilibration.inp step5_production.inp | tr '\n' ' ')"

echo "System preparation complete."
exit 0
