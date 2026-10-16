#!/usr/bin/env bash
# deploy_nira_prot.sh — niraparib protonated rerun (S1 24 ns + S2 26 ns)
# Responds to reviewer comment R3-M8: the piperidine N (pKa 9.95) is protonated (+1) at pH 7.4; the original model was neutral.
# Usage: on the AutoDL instance (conda env `gamd`): bash md/autodl_setup/deploy_nira_prot.sh
# Prerequisite: already uploaded (see upload_nira_prot.sh):
#   docking/results_system2/mol2/niraparib_best_prot.mol2   (S2 pose, protonated, rdkit-converted)
#   docking/results/niraparib_docked_prot.sdf               (S1 pose, protonated)
#   data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb (S2 receptor)
#   docking/receptors/6VKK_A_receptor.pdb                   (S1 receptor)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REC_S1="$ROOT/docking/receptors/6VKK_A_receptor.pdb"
LIG_S1="$ROOT/docking/results/niraparib_docked_prot.sdf"
REC_S2="$ROOT/data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb"
LIG_S2="$ROOT/docking/results_system2/mol2/niraparib_best_prot.mol2"
RUNS="$ROOT/runs_nira_prot"
ARCHIVE="$ROOT/results/rerun_logs"
mkdir -p "$RUNS" "$ARCHIVE"

PROD_S2=26
PROD_S1=24
SEED=42
TEMP=300
NIRA_PROT_TAG="sys2_niraparib_prot"

# ---------- 1. S2: tleap build (internal antechamber -nc 1, patched) ----------
echo "==== [1/6] S2 tleap build (protonated niraparib) ===="
python "$HERE/build_system2_tleap.py" \
    --receptor "$REC_S2" --ligand "$LIG_S2" --out "$RUNS/$NIRA_PROT_TAG" 2>&1 | tail -5

# ---------- 3. S2: equilibration + GaMD 26 ns ----------
echo "==== [3/6] S2 equilibrate + GaMD ${PROD_S2} ns ===="
S2PRM="$RUNS/$NIRA_PROT_TAG/${NIRA_PROT_TAG}.prmtop"
S2INC="$RUNS/$NIRA_PROT_TAG/${NIRA_PROT_TAG}.inpcrd"
CUDA_VISIBLE_DEVICES=0 python "$HERE/equilibrate.py" --prmtop "$S2PRM" --inpcrd "$S2INC" \
    --out "$RUNS/$NIRA_PROT_TAG/${NIRA_PROT_TAG}_equil" --temp $TEMP --nvt-ns 0.2 --npt-ns 1.0 2>&1 | tail -3
python "$HERE/gen_gamd_config.py" --prmtop "$S2PRM" \
    --rst7 "$RUNS/$NIRA_PROT_TAG/${NIRA_PROT_TAG}_equil.rst7" \
    --outdir "$RUNS/$NIRA_PROT_TAG/gamd_out" --out "$RUNS/$NIRA_PROT_TAG/config.xml" \
    --prod-ns $PROD_S2 --no-min --boost-type lower-dual --seed $SEED \
    --report-ps 50 --temp $TEMP 2>&1 | tail -2
CUDA_VISIBLE_DEVICES=0 nohup "$ROOT/tools/gamd-openmm/gamdRunner" xml \
    "$RUNS/$NIRA_PROT_TAG/config.xml" > "$RUNS/$NIRA_PROT_TAG/runner.log" 2>&1 < /dev/null &
echo "  [S2] gamdRunner PID $! -> GPU 0"

# ---------- 4. S1: build (OpenFF gaff-2.11 auto-parameterisation; protonated sdf keeps +1) ----------
echo "==== [4/6] S1 build (OpenFF, protonated sdf) ===="

# ---------- 5. S1: build + equilibration + GaMD 24 ns ----------
echo "==== [5/6] S1 build + equilibrate + GaMD ${PROD_S1} ns ===="
S1TAG="sys1_niraparib_prot"
python "$HERE/build_system.py" --receptor "$REC_S1" --ligand "$LIG_S1" --out "$RUNS/$S1TAG" 2>&1 | tail -4
S1PRM="$RUNS/$S1TAG.prmtop"
S1INC="$RUNS/$S1TAG.inpcrd"
CUDA_VISIBLE_DEVICES=1 python "$HERE/equilibrate.py" --prmtop "$S1PRM" --inpcrd "$S1INC" \
    --out "$RUNS/$S1TAG/${S1TAG}_equil" --temp $TEMP --nvt-ns 0.2 --npt-ns 1.0 2>&1 | tail -3
python "$HERE/gen_gamd_config.py" --prmtop "$S1PRM" \
    --rst7 "$RUNS/$S1TAG/${S1TAG}_equil.rst7" \
    --outdir "$RUNS/$S1TAG/gamd_out" --out "$RUNS/$S1TAG/config.xml" \
    --prod-ns $PROD_S1 --no-min --boost-type lower-dual --seed $SEED \
    --report-ps 50 --temp $TEMP 2>&1 | tail -2
CUDA_VISIBLE_DEVICES=1 nohup "$ROOT/tools/gamd-openmm/gamdRunner" xml \
    "$RUNS/$S1TAG/config.xml" > "$RUNS/$S1TAG/runner.log" 2>&1 < /dev/null &
echo "  [S1] gamdRunner PID $! -> GPU 1"

echo "==== [6/6] both systems launched (S2 26 ns GPU0, S1 24 ns GPU1) ===="
echo "monitor: tail -f $RUNS/$NIRA_PROT_TAG/runner.log | tail -3"
echo "estimated S2 ~17h (285K atoms ~70 ns/day), S1 ~2h (61K atoms)"
