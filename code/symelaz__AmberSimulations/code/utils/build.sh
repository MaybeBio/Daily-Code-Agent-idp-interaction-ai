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
    a99SBdisp_tip4pd)
        WATER=${WATER:-authors}
        case "$WATER" in
            authors) WATERBOX=TIP4PEWBOX ;;
            amber)   WATERBOX=TIP4PDBOX ;;
            *) echo "ERROR: WATER must be 'authors' or 'amber'"; exit 1 ;;
        esac
        LEAPRC=leaprc.pipe_a99SBdisp_tip4pd_${WATER}
        SALT=${SALT:-0.15}; CATION=${CATION:-K+};   ADD8=true ;;
    *) echo "ERROR: MODEL must be 'ff19SB_tip3p' or 'a99SBdisp_tip4pd'"; exit 1 ;;
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
        -e "s/__CATION__/$CATION/g" -e "s/__SALT__/$SALT/" -e "s/__NSALT__/$NSALT/g" -e "s/__SALTLINE__/$SALTLINE/" "$1"
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
    antechamber -i "$FOLDER"/out.mol2 -fi mol2 -o "$FOLDER"/out.pdb -fo pdb

    LEAPFILE=utils/build_complex.leap
else
    echo "No ligand detected. Only protein will be processed."
    LEAPFILE=utils/build_protein.leap
fi

# Extract protein without hydrogens
vmd -dispdev text -e utils/extract_protein.tcl -args "$PDBFILE" "$FOLDER/protein_noh.pdb" > $FOLDER/extract_protein.log

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
python utils/check_topology.py step3_input.parm7 $MODEL $WATER | tee $FOLDER/topology_report.txt
if grep -q "^FAIL" $FOLDER/topology_report.txt; then
    echo "ERROR: topology check failed (see $FOLDER/topology_report.txt)"; exit 1
fi

# OpenMM inputs for this model (edit these copies, not the templates in inputs/)
rm -f step4.*_equilibration.inp step5_production.inp
cp inputs/$MODEL/*.inp .
echo "Copied OpenMM inputs from inputs/$MODEL/: $(ls step4.*_equilibration.inp step5_production.inp | tr '\n' ' ')"

echo "System preparation complete."
exit 0
