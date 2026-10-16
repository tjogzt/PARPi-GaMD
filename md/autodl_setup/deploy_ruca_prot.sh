#!/usr/bin/env bash
# deploy_ruca_prot.sh — rucaparib protonated rerun (S1 24 ns + S2 26 ns)
# Responds to external review 3.1/R3: rucaparib benzylamine (CH2NHCH3) pKa ~8.4, ~9% protonated (+1) at pH 7.4.
# Same protocol as the niraparib protonation: lower-dual, sigma0 6.0, dt 2 fs, seed 42, 300 K,
#   S1 24 ns / S2 26 ns, report-ps 50, production mask 6.0e6 steps, exact DBE reweighting.
# Usage: on the AutoDL instance (conda env `gamd`): bash md/autodl_setup/deploy_ruca_prot.sh
# Prerequisite: already uploaded (upload_ruca_prot.sh):
#   docking/results_system2/mol2/rucaparib_best_prot.mol2   (S2 pose, protonated, +1)
#   docking/results/rucaparib_docked_prot.sdf               (S1 pose, protonated, +1)
#   data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb (S2 receptor)
#   docking/receptors/6VKK_A_receptor.pdb                   (S1 receptor)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REC_S1="$ROOT/docking/receptors/6VKK_A_receptor.pdb"
LIG_S1="$ROOT/docking/results/rucaparib_docked_prot.sdf"
REC_S2="$ROOT/data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb"
LIG_S2="$ROOT/docking/results_system2/mol2/rucaparib_best_prot.mol2"
RUNS="$ROOT/runs_ruca_prot"
ARCHIVE="$ROOT/results/rerun_logs"
mkdir -p "$RUNS" "$ARCHIVE"

PROD_S2=26
PROD_S1=24
SEED=42
TEMP=300
S2TAG="sys2_rucaparib_prot"
S1TAG="sys1_rucaparib_prot"

# ---------- 1. S2: tleap build (antechamber -nc 1; outputs under <out>/build/) ----------
echo "==== [1/6] S2 tleap build (protonated rucaparib) ===="
python "$HERE/build_system2_tleap.py" \
    --receptor "$REC_S2" --ligand "$LIG_S2" --out "$RUNS/$S2TAG" 2>&1 | tail -4

# ---------- 2. S2: equilibration + GaMD 26 ns (GPU 0) ----------
echo "==== [2/6] S2 equilibrate + GaMD ${PROD_S2} ns ===="
S2PRM="$RUNS/$S2TAG/build/rucaparib_best_prot.prmtop"
S2INC="$RUNS/$S2TAG/build/rucaparib_best_prot.inpcrd"
test -f "$S2PRM" || { echo "FATAL: $S2PRM missing"; exit 1; }
CUDA_VISIBLE_DEVICES=0 python "$HERE/equilibrate.py" --prmtop "$S2PRM" --inpcrd "$S2INC" \
    --out "$RUNS/$S2TAG/${S2TAG}_equil" --temp $TEMP --nvt-ns 0.2 --npt-ns 1.0 2>&1 | tail -3
python "$HERE/gen_gamd_config.py" --prmtop "$S2PRM" \
    --rst7 "$RUNS/$S2TAG/${S2TAG}_equil.rst7" \
    --outdir "$RUNS/$S2TAG/gamd_out" --out "$RUNS/$S2TAG/config.xml" \
    --prod-ns $PROD_S2 --no-min --boost-type lower-dual --seed $SEED \
    --report-ps 50 --temp $TEMP 2>&1 | tail -2
CUDA_VISIBLE_DEVICES=0 nohup "$ROOT/tools/gamd-openmm/gamdRunner" xml \
    "$RUNS/$S2TAG/config.xml" > "$RUNS/$S2TAG/runner.log" 2>&1 < /dev/null &
echo "  [S2] gamdRunner PID $! -> GPU 0"

# ---------- 3. S1: build (OpenFF gaff-2.11; protonated sdf keeps +1) ----------
echo "==== [3/6] S1 build (OpenFF, protonated sdf) ===="
python "$HERE/build_system.py" --receptor "$REC_S1" --ligand "$LIG_S1" --out "$RUNS/$S1TAG" 2>&1 | tail -4
S1PRM="$RUNS/$S1TAG.prmtop"
S1INC="$RUNS/$S1TAG.inpcrd"
test -f "$S1PRM" || { echo "FATAL: $S1PRM missing"; exit 1; }

# ---------- 4. S1: equilibration + GaMD 24 ns (GPU 1) ----------
echo "==== [4/6] S1 equilibrate + GaMD ${PROD_S1} ns ===="
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

# ---------- 5. done ----------
echo "==== [5/5] both systems launched (S2 26 ns GPU0, S1 24 ns GPU1) ===="
echo "monitor: tail -f $RUNS/$S2TAG/runner.log | tail -3"
echo "estimated S2 ~17h (285K atoms ~70 ns/day), S1 ~2h (61K atoms)"
