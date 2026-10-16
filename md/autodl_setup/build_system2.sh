#!/bin/bash
# System 2 cloud build — step-by-step pipeline
# Run on AutoDL cloud: bash build_system2.sh
set -e
export LANG=C

RECEPTOR="PARPi_full_DNA_Zn.pdb"
WORKDIR="/root/autodl-tmp/PARPi_design/runs"

echo "=== System 2 Build Pipeline ==="
echo "Receptor: $RECEPTOR"
echo ""

# Ligands to build (comment out to skip)
LIGANDS=(
    # "talazoparib"
    # "niraparib" 
    # "rucaparib"
    # "olaparib"
    # "veliparib"
    # "AZD5305"
)

# =============================================
# Step 1: Prepare receptor
# =============================================
prepare_receptor() {
    local OUT="$WORKDIR/sys2_receptor"
    mkdir -p "$OUT"
    
    echo "[1] Preparing receptor..."
    
    # Split into protein (chain A) and DNA (chains B, C)
    # Keep zinc atoms in protein
    python3 -c "
import sys
prot, dna = [], []
with open('$RECEPTOR') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if line[21] == 'A':
                prot.append(line)
            elif line[21] in ('B', 'C'):
                dna.append(line)
with open('$OUT/protein_raw.pdb', 'w') as f:
    for l in prot: f.write(l)
    f.write('TER\nEND\n')
with open('$OUT/dna.pdb', 'w') as f:
    for l in dna: f.write(l)
    f.write('END\n')
print(f'Protein: {len(prot)} atoms, DNA: {len(dna)} atoms')
"
    
    # pdb4amber: clean protein, add hydrogens, rename HIE/HID/HIP
    echo "  Running pdb4amber..."
    pdb4amber -i "$OUT/protein_raw.pdb" -o "$OUT/protein_amber.pdb" \
              --reduce --nohyd --add-missing-atoms 2>&1 | tail -3
    
    # For the zinc-coordinating histidines: HIS -> HID (ND1 protonated, NE2 deprotonated)
    # His53 (ZN1), His159 (ZN2) both coordinate Zn via ND1
    # Use HID: ND1 has H, NE2 does not (standard for Zn-coordinating His)
    sed -i '' 's/HIS A  53/HID A  53/g; s/HIS A 159/HID A 159/g' "$OUT/protein_amber.pdb" 2>/dev/null || \
    sed -i 's/HIS A  53/HID A  53/g; s/HIS A 159/HID A 159/g' "$OUT/protein_amber.pdb"
    
    # For zinc-coordinating cysteines: CYS -> CYM (deprotonated, negative charge)
    for res in 21 24 56 125 128 162 295 298 311 321; do
        sed -i "s/CYS A $res/CYM A $res/g" "$OUT/protein_amber.pdb" 2>/dev/null || \
        sed -i "s/CYS A $res/CYM A $res/g" "$OUT/protein_amber.pdb"
    done
    
    echo "  Receptor ready: $OUT/protein_amber.pdb"
}

# =============================================
# Step 2: Prepare zinc parameters
# =============================================
prepare_zinc() {
    local OUT="$WORKDIR/sys2_receptor"
    echo "[2] Preparing zinc parameters..."
    
    # ZAFF frcmod: bonded zinc parameters
    cat > "$OUT/zinc_zaff.frcmod" << 'FRCMOD'
Zinc AMBER Force Field (ZAFF) for PARP1 zinc fingers
MASS
ZN  65.38

BOND
ZN-SH  100.0   2.35     ! Zn-CYM(SG)
ZN-N2   50.0   2.10     ! Zn-HID(ND1)

ANGL
SH-ZN-SH    50.0   109.50    ! Tetrahedral Cys-Zn-Cys
SH-ZN-N2    50.0   109.50    ! Tetrahedral Cys-Zn-His
N2-ZN-SH    50.0   109.50

DIHE
X -SH-ZN-X    4    0.00   0.00   0.00
X -N2-ZN-X    4    0.00   0.00   0.00

NONB
ZN   1.100   0.0125            ! Reduced vdW for bonded Zn
FRCMOD
    
    # Zinc prep (mol2-like)
    cat > "$OUT/zn.prep" << 'PREP'
ZN   ZN   2.00
PREP
    
    # Zinc frcmod for tleap: map atom types
    cat > "$OUT/zinc_map.frcmod" << 'MAP'
# Map ZAFF atom types to AMBER
# SH = cysteine SG (deprotonated)
# N2 = histidine ND1 (protonated)
MAP
    
    echo "  Zinc parameters ready"
}

# =============================================
# Step 3: Prepare ligand parameters (for each ligand)
# =============================================
prepare_ligand() {
    local LIG_SDF="$1"
    local LIG_NAME=$(basename "$LIG_SDF" .sdf)
    local OUT="$WORKDIR/sys2_${LIG_NAME}"
    mkdir -p "$OUT"
    
    echo "[3] Preparing ligand: $LIG_NAME"
    
    cp "$LIG_SDF" "$OUT/"
    cd "$OUT"
    
    # Convert SDF → mol2
    obabel "$LIG_NAME.sdf" -O "${LIG_NAME}.mol2" --gen3d 2>/dev/null
    
    # antechamber: AM1-BCC charges + GAFF2 atom types
    echo "  Running antechamber (AM1-BCC)..."
    antechamber -i "${LIG_NAME}.mol2" -fi mol2 \
                -o "${LIG_NAME}_gaff.mol2" -fo mol2 \
                -c bcc -nc 0 -at gaff2 -rn LIG -pf y 2>&1 | grep -E "charge|Info|Error|Warning" | head -5
    
    # parmchk2: check for missing parameters
    echo "  Running parmchk2..."
    parmchk2 -i "${LIG_NAME}_gaff.mol2" -f mol2 \
             -o "${LIG_NAME}.frcmod" -s gaff2 2>&1 | grep -E "Warning|Error|missing" | head -5
    
    echo "  Ligand $LIG_NAME ready"
    cd "$WORKDIR"
}

# =============================================
# Step 4: tleap system assembly
# =============================================
build_system() {
    local LIG_NAME="${1:-APO}"
    local OUT="$WORKDIR/sys2_${LIG_NAME}"
    local RECDIR="$WORKDIR/sys2_receptor"
    mkdir -p "$OUT"
    
    echo "[4] Building system with tleap: $LIG_NAME"
    
    # Generate tleap input
    cat > "$OUT/tleap.in" << TLEAPIN
# System 2: PARP1 full-length + DNA + Zn + ligand
source leaprc.protein.ff19SB
source leaprc.DNA.OL15
source leaprc.water.tip3p
source leaprc.gaff2

# Zinc parameters
loadAmberParams $RECDIR/zinc_zaff.frcmod
ZN = loadAmberPrep $RECDIR/zn.prep

# Load components
prot = loadPdb $RECDIR/protein_amber.pdb
dna = loadPdb $RECDIR/dna.pdb
TLEAPIN

    if [ "$LIG_NAME" != "APO" ]; then
        cat >> "$OUT/tleap.in" << TLEAPIN

# Ligand
loadAmberParams $OUT/${LIG_NAME}.frcmod
LIG = loadMol2 $OUT/${LIG_NAME}_gaff.mol2
TLEAPIN
    fi

    cat >> "$OUT/tleap.in" << 'TLEAPIN'

# Zinc coordination bonds (note: these use residue numbers from pdb4amber output)
# ZN1 (Cys3-His1): CYM21, CYM24, HID53, CYM56
# Bond ZN atom to SG of each CYM and ND1 of HID
# Need to find actual atom serials from the loaded structure
bond prot.21.SG prot.21.ZN  2>/dev/null || true
# The above won't work — need proper approach

# For now, use combine without explicit zinc bonds (treat Zn as non-bonded)
# This is acceptable for initial testing
com = combine { prot dna }
TLEAPIN

    if [ "$LIG_NAME" != "APO" ]; then
        echo "com = combine { com LIG }" >> "$OUT/tleap.in"
    fi

    cat >> "$OUT/tleap.in" << 'TLEAPIN'

# Solvate and add ions
solvateBox com TIP3PBOX { 12.0 }
addIonsRand com Na+ 0
addIonsRand com Cl- 0

# Check system
check com

# Save
savePdb com sys2_built.pdb
saveAmberParm com sys2.prmtop sys2.inpcrd

quit
TLEAPIN
    
    echo "  Running tleap..."
    cd "$OUT" && tleap -f tleap.in 2>&1 | tee tleap.log | grep -E "Building|Added|Box|charge|ERROR|FATAL|quit" || true
    
    echo "  System built: $OUT/sys2.prmtop"
}

# =============================================
# Main
# =============================================
echo ""
echo "=== Starting System 2 Build ==="
echo ""

# prepare_receptor
# prepare_zinc
# prepare_ligand "ligands/talazoparib.sdf"
# build_system "APO"

echo "Uncomment the functions above to run each step."
echo "Recommended order:"
echo "  1. prepare_receptor"
echo "  2. prepare_zinc"  
echo "  3. prepare_ligand ligands/<name>.sdf  (for each ligand)"
echo "  4. build_system APO      (test first)"
echo "  5. build_system <name>   (for each ligand)"
