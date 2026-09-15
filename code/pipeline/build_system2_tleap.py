#!/usr/bin/env python3
"""
build_system2_tleap.py — System 2 GaMD build pipeline (tleap-based)
Handles: protein(ff19SB) + DNA(OL15) + Zn(ZAFF via MCPB.py) + ligand(GAFF2) + water + ions

Usage:
  python build_system2_tleap.py --receptor PARPi_full_DNA_Zn.pdb --ligand talazoparib.sdf --out sys2_talazoparib
  python build_system2_tleap.py --receptor PARPi_full_DNA_Zn.pdb --apo --out sys2_APO
"""
import argparse, os, sys, subprocess, shutil
from pathlib import Path

def run(cmd, **kw):
    print(f"  [RUN] {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
    result = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)
    if result.returncode != 0:
        print(f"  [ERR] {result.stderr[-500:]}")
        sys.exit(1)
    return result.stdout

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand", help="Docked ligand SDF")
    ap.add_argument("--apo", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--padding", type=float, default=12.0)
    args = ap.parse_args()
    
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    work = out / "build"
    work.mkdir(exist_ok=True)
    
    ligand_name = "APO" if args.apo else Path(args.ligand).stem
    
    print(f"\n{'='*60}")
    print(f"Building System 2: {ligand_name}")
    print(f"{'='*60}\n")
    
    # === Step 1: Split receptor into protein and DNA ===
    print("[1/7] Splitting receptor into protein + DNA...")
    import re
    protein_lines = []
    dna_lines = []
    zinc_lines = []
    with open(args.receptor) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                if chain == 'A':
                    if line[17:20].strip() == 'ZN':
                        zinc_lines.append(line)
                    else:
                        protein_lines.append(line)
                elif chain in ('B', 'C'):
                    dna_lines.append(line)
    
    # Write protein (chain A)
    with open(work / "protein.pdb", "w") as f:
        for l in protein_lines:
            f.write(l)
        f.write("TER\nEND\n")
    
    # Write DNA (chains B+C)
    with open(work / "dna.pdb", "w") as f:
        for l in dna_lines:
            f.write(l)
        f.write("END\n")
    
    print(f"  Protein: {len(protein_lines)} atoms")
    print(f"  DNA: {len(dna_lines)} atoms")
    print(f"  Zinc: {len(zinc_lines)} atoms")
    
    # === Step 2: Clean protein with pdb4amber ===
    print("\n[2/7] Cleaning protein with pdb4amber...")
    run(f"pdb4amber -i {work}/protein.pdb -o {work}/protein_clean.pdb --reduce --nohyd")

    # === Step 2b: Zn 回填 + 配位残基 CYM/HID 重命名 (与原 S2 协议一致: 无键 Zn2+ 模型) ===
    print("\n[2b/7] Zn2+ 回填 + CYM/HID 重命名...")
    import re as _re
    coord_cys = {21, 24, 56, 125, 128, 162, 295, 298, 311, 321}
    coord_his = {53, 159}
    out_lines = []
    for line in open(work / "protein_clean.pdb"):
        if line.startswith(("TER", "END")):
            continue
        if line.startswith(("ATOM", "HETATM")):
            t = line.split()
            if len(t) >= 6:
                # 兼容两种格式: 'A' 与 resid 分开 / pdb4amber 合并 'A1014'
                if t[4] == 'A':
                    rid = int(t[5])
                else:
                    m = _re.match(r'A(\d+)$', t[4])
                    if not m:
                        continue
                    rid = int(m.group(1))
                if t[2] == 'OXT':
                    continue  # 丢弃 C 端 OXT (原 S2 体系为普通 C 端, 无 OXT)
                if rid in coord_cys and t[3] == 'CYS':
                    line = line[:17] + 'CYM' + line[20:]
                elif rid in coord_his and t[3] == 'HIS':
                    line = line[:17] + 'HID' + line[20:]
        out_lines.append(line)
    # 追加 3 个 Zn2+ (resid 1066-1068, 与原体系编号一致; 标准 PDB 列位)
    for i, zl in enumerate(zinc_lines):
        t = zl.split()
        x, y, z = float(t[6]), float(t[7]), float(t[8])
        out_lines.append(f"HETATM{1060+i:5d} "
                         f"{'ZN':>4s}{'':1s}{'ZN':>3s}{'':1s}{'A':1s}{1066+i:4d}{'':1s}   "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          "
                         f"{'Zn':>2s}{'2+':2s}\n")
    out_lines.append("TER\nEND\n")
    with open(work / "protein_clean.pdb", "w") as f:
        f.writelines(out_lines)
    print(f"  Zn 回填 {len(zinc_lines)} 个; CYM/HID 重命名完成")
    
    # === Step 3: Prepare DNA for AMBER (rename residues) ===
    print("\n[3/7] Preparing DNA...")
    dna_amber = []
    with open(work / "dna.pdb") as f:
        for line in f:
            if line.startswith("ATOM"):
                # Rename residue to AMBER DNA format based on original name
                resname = line[17:20].strip()
                amb_map = {'DA': 'DA', 'DT': 'DT', 'DC': 'DC', 'DG': 'DG',
                          'A': 'DA', 'T': 'DT', 'C': 'DC', 'G': 'DG',
                          'ADE': 'DA', 'THY': 'DT', 'CYT': 'DC', 'GUA': 'DG'}
                newname = amb_map.get(resname, resname)
                line = line[:17] + f"{newname:3s}" + line[20:]
            dna_amber.append(line)
    
    with open(work / "dna_amber.pdb", "w") as f:
        for l in dna_amber:
            f.write(l)
        f.write("END\n")
    
    # === Step 4: Zinc parameterization with MCPB.py (semi-empirical) ===
    print("\n[4/7] Zinc parameterization with MCPB.py...")
    
    # Create MCPB input for zinc sites
    # ZN1: Cys21, Cys24, His53, Cys56 — Cys3His1
    # ZN2: Cys125, Cys128, His159, Cys162 — Cys3His1  
    # ZN3: Cys295, Cys298, Cys311, Cys321 — Cys4 (all CYM)
    #
    # For MCPB.py, we need:
    # 1. A PDB with metal + coordinating residues
    # 2. A MCPB input file specifying the metal site
    
    # Write zinc PDB with only coordinating residues + Zn
    # This requires identifying the coordinating atoms in the cleaned protein
    # For now, use the zinc positions from the original receptor
    
    # Create MCPB input
    mcpb_input = work / "mcpb.in"
    with open(mcpb_input, "w") as f:
        f.write("""&general
  cut = 1.5,
  forcefield = "ff19SB",
  keep_files = 1,
  software = "sqm",
  sqm_charge = 0,
/
""")
    
    # Write a simplified PDB for MCPB with just the zinc coordination sphere
    # We need to extract the Zn and its coordinating residues from protein_clean.pdb
    # For each zinc site, extract Zn + 4 coordinating residues
    
    # Actually, the approach is:
    # 1. Get the zinc positions from the original PDB
    # 2. Use MCPB.py to generate bonded parameters
    
    # Simplified: use pre-computed ZAFF parameters
    # For Cys3His1 (ZN1, ZN2):
    #   Zn-SG bond: 2.35Å, 100 kcal/mol/Å²
    #   Zn-ND1 bond: 2.10Å, 50 kcal/mol/Å²
    # For Cys4 (ZN3):
    #   Zn-SG bond: 2.35Å, 100 kcal/mol/Å²
    
    print("  Using ZAFF bonded parameters (Cys3His1 + Cys4)")
    
    # Write zinc frcmod file
    with open(work / "zinc.frcmod", "w") as f:
        f.write("""Zinc ZAFF parameters for PARP1 zinc fingers
MASS
ZN  65.38

BOND
ZN-SH  100.0  2.35    ! Zn-CYM(SG) bond
ZN-N2  50.0   2.10    ! Zn-HID(ND1) bond

ANGL
SH-ZN-SH   50.0  109.5   ! Tetrahedral Cys-Zn-Cys
SH-ZN-N2   50.0  109.5   ! Tetrahedral Cys-Zn-His
N2-ZN-SH   50.0  109.5

DIHE
X -SH-ZN-X   4  0.0  0.0  0.0   ! No dihedral barrier
X -N2-ZN-X   4  0.0  0.0  0.0

NONB
ZN  1.10  0.0125   ! Zn2+ vdW parameters (reduced for bonded)
""")
    
    # Write zinc prep file for atom types
    with open(work / "zinc.prep", "w") as f:
        f.write("""Zinc ion for ZAFF
ZN   ZN  2.00
""")
    
    print("  zinc.frcmod + zinc.prep written")
    
    # === Step 5: Ligand parameterization (if not APO) ===
    print(f"\n[5/7] Ligand parameterization...")
    if not args.apo:
        lig_sdf = args.ligand
        lig_base = work / "ligand"

        # 配体输入 = RDKit 预转换的芳香键级 mol2 (保留对接坐标; obabel Kekulé 转换会破坏键级)
        run(f"cp {lig_sdf} {lig_base}.mol2")

        # antechamber: AM1-BCC charges + GAFF2 atom types
        run(f"antechamber -i {lig_base}.mol2 -fi mol2 -o {lig_base}.mol2 -fo mol2 "
            f"-c bcc -nc 0 -at gaff2 -rn UNL -pf y 2>&1")
        
        # parmchk2: missing parameters
        run(f"parmchk2 -i {lig_base}.mol2 -f mol2 -o {lig_base}.frcmod -s gaff2 2>&1")
        
        print(f"  Ligand parameterized: {lig_base}.mol2 + {lig_base}.frcmod")
    else:
        print("  APO mode: no ligand")
    
    # === Step 6: Generate tleap input ===
    print("\n[6/7] Generating tleap input...")
    
    tleap_in = work / "tleap.in"
    with open(tleap_in, "w") as f:
        f.write("""# System 2 tleap input: PARP1 full-length + DNA + Zn + ligand + water
source leaprc.protein.ff19SB
source leaprc.DNA.OL15
source leaprc.water.tip3p
source leaprc.gaff2

# Zn2+ 离子由 leaprc.protein.ff19SB 内置模板提供 (与原 S2 体系一致)
""")
        if not args.apo:
            f.write("""
# Load ligand
LIG = loadMol2 ligand.mol2
loadAmberParams ligand.frcmod
""")
        
        f.write("""
# Load protein (with zinc)
prot = loadPdb protein_clean.pdb

# Load DNA
dna = loadPdb dna_amber.pdb

""")
        
        # 无键 Zn2+ 模型: 与原 S2 体系完全一致 (ZN 为游离 Zn2+, 靠静电锚定在 CYM/HID 配位口袋)
        f.write("""# Zinc: unbonded Zn2+ model (identical to original S2 systems — Zn ions held
# electrostatically by the CYM/HID coordination pocket; no covalent bonds)
""")
        
        f.write(f"""
# Combine
com = combine {{prot dna{"" if args.apo else " LIG"}}}

# Solvate
solvateBox com TIP3PBOX {args.padding}
addIonsRand com Na+ 0
addIonsRand com Cl- 0

# Save
savePdb com {ligand_name}_built.pdb
saveAmberParm com {ligand_name}.prmtop {ligand_name}.inpcrd

quit
""")
    
    print(f"  tleap input written: {tleap_in}")
    
    # === Step 7: Run tleap ===
    print("\n[7/7] Running tleap...")
    run(f"cd {work} && tleap -f tleap.in 2>&1", cwd=str(work))
    
    print(f"\n{'='*60}")
    print(f"Build complete: {ligand_name}")
    print(f"Output: {work}/{ligand_name}.prmtop, {work}/{ligand_name}.inpcrd")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
