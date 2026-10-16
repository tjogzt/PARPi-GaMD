#!/usr/bin/env bash
# deploy_veli_prot_resume.sh — resume after interruption: the S2 tleap build already finished (21:39); start from equilibration.
# Usage: on the instance (conda env `gamd`, GPU required): bash md/autodl_setup/deploy_veli_prot_resume.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REC_S1="$ROOT/docking/receptors/6VKK_A_receptor.pdb"
LIG_S1="$ROOT/docking/results/veliparib_docked_prot.sdf"
RUNS="$ROOT/runs_veli_prot"

PROD_S2=26
PROD_S1=24
SEED=42
TEMP=300
S2TAG="sys2_veliparib_prot"
S1TAG="sys1_veliparib_prot"

# ---------- 2. S2: equilibration + GaMD 26 ns (GPU 0) ----------
echo "==== [2/6-resume] S2 equilibrate + GaMD ${PROD_S2} ns ===="
S2PRM="$RUNS/$S2TAG/build/veliparib_best_prot.prmtop"
S2INC="$RUNS/$S2TAG/build/veliparib_best_prot.inpcrd"
test -f "$S2PRM" || { echo "FATAL: $S2PRM missing (build not complete)"; exit 1; }
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
echo "==== [3/6-resume] S1 build (OpenFF, protonated sdf) ===="
python "$HERE/build_system.py" --receptor "$REC_S1" --ligand "$LIG_S1" --out "$RUNS/$S1TAG" 2>&1 | tail -4
S1PRM="$RUNS/$S1TAG.prmtop"
S1INC="$RUNS/$S1TAG.inpcrd"
test -f "$S1PRM" || { echo "FATAL: $S1PRM missing"; exit 1; }

# ---------- 4. S1: equilibration + GaMD 24 ns (GPU 1) ----------
echo "==== [4/6-resume] S1 equilibrate + GaMD ${PROD_S1} ns ===="
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

echo "==== [5/5] both systems launched (S2 26 ns GPU0, S1 24 ns GPU1) ===="
