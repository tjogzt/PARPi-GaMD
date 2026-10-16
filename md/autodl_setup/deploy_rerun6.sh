#!/usr/bin/env bash
# deploy_rerun6.sh — rerun 6 systems to replace reconstructed boost energies
# with exact archived logs (removes the "rec." flags in Table S5 and the S1
# rerun backlog from the 2026-09-21 reweighting audit).
#
#   Batch 1 (S2, three original-panel drugs, 26 ns production):
#       reuse the archived prmtop + the frame-0 restart extracted from the
#       original production DCD (md/results/rerun_starts/sys2_<lig>_start.rst7)
#       -> gamdRunner re-equilibrates (cMD 2 ns + GaMD eq 10 ns) then produces
#       26 ns. gamd.log IS archived this time (see ARCHIVE dir).
#   Batch 2 (S1, three original-panel drugs, 24 ns production):
#       full rebuild from scratch (CAT-domain receptor 6VKK_A + docked pose)
#       -> equilibrate -> GaMD 24 ns.
#
# Layout: runs_rerun/sys1_<lig>/, runs_rerun/sys2_<lig>/
# Usage:   bash deploy_rerun6.sh    (on the AutoDL instance, conda env `gamd`)
# Requires: 3 GPUs (CUDA_VISIBLE_DEVICES 0/1/2), ~1 day total.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REC_S1="$ROOT/docking/receptors/6VKK_A_receptor.pdb"
LIG_S1="$ROOT/docking/results"
REC_S2_DATA="/Volumes/tjogzt4T/PARPi_data"   # source of archived S2 prmtops
RUNS="$ROOT/runs_rerun"
ARCHIVE="$ROOT/results/rerun_logs"
mkdir -p "$RUNS" "$ARCHIVE"

# On the instance, the archived S2 prmtops and extracted restart files are
# uploaded to md/rerun_inputs/ (see upload script). If the local data volume
# is visible (dev machine), fall back to it.
S2_INPUTS="$ROOT/md/rerun_inputs"
[ -d "$REC_S2_DATA" ] && [ ! -d "$S2_INPUTS" ] && S2_INPUTS="$REC_S2_DATA"

PROD_S2=26
PROD_S1=24
SEED=42
TEMP=300

echo "==== [S2] Batch 1: niraparib/olaparib/rucaparib (26 ns, reused topologies) ===="
for drug in niraparib olaparib rucaparib; do
  TAG="sys2_${drug}"
  PRM="$RUNS/$TAG/sys2_${drug}.prmtop"
  RST="$RUNS/$TAG/sys2_${drug}_start.rst7"
  mkdir -p "$RUNS/$TAG"
  cp "$S2_INPUTS/sys2_${drug}/sys2_${drug}.prmtop" "$PRM"
  cp "$S2_INPUTS/sys2_${drug}/sys2_${drug}_start.rst7" "$RST"
  python "$HERE/gen_gamd_config.py" --prmtop "$PRM" --rst7 "$RST" \
      --outdir "$RUNS/$TAG/gamd_out" --out "$RUNS/$TAG/config.xml" \
      --prod-ns $PROD_S2 --no-min --boost-type lower-dual --seed $SEED \
      --report-ps 50 --temp $TEMP 2>&1 | tail -2
  echo "  [$TAG] config ready"
done

# launch S2 on GPUs 0/1/2
i=0
for drug in niraparib olaparib rucaparib; do
  TAG="sys2_${drug}"
  rm -f "$RUNS/$TAG/gamd_out/gamd.log" "$RUNS/$TAG/gamd_out/output.dcd"
  CUDA_VISIBLE_DEVICES=$i nohup "$ROOT/tools/gamd-openmm/gamdRunner" xml \
      "$RUNS/$TAG/config.xml" > "$RUNS/$TAG/runner.log" 2>&1 < /dev/null &
  echo "  [$TAG] gamdRunner PID $! -> GPU $i"
  i=$((i+1))
  sleep 5
done
echo "==== S2 batch launched; waiting for completion ===="
# wait for all three runner.log to contain the completion marker or processes to end
while pgrep -f gamdRunner > /dev/null; do sleep 600; done
for drug in niraparib olaparib rucaparib; do
  TAG="sys2_${drug}"
  cp "$RUNS/$TAG/gamd_out/gamd.log" "$ARCHIVE/sys2_${drug}_gamd.log"
  echo "  [$TAG] gamd.log archived ($(wc -l < "$ARCHIVE/sys2_${drug}_gamd.log") lines)"
done

echo "==== [S1] Batch 2: niraparib/olaparib/rucaparib (24 ns, full rebuild) ===="
for drug in niraparib olaparib rucaparib; do
  TAG="sys1_${drug}"
  PRM="$RUNS/${TAG}.prmtop"
  INC="$RUNS/${TAG}.inpcrd"
  if [ ! -f "$PRM" ]; then
    python "$HERE/build_system.py" --receptor "$REC_S1" \
        --ligand "$LIG_S1/${drug}_docked.fixed.sdf" --out "$RUNS/$TAG" 2>&1 | tail -4
  fi
  CUDA_VISIBLE_DEVICES=2 python "$HERE/equilibrate.py" \
      --prmtop "$PRM" --inpcrd "$INC" --out "$RUNS/$TAG/sys1_${drug}_equil" \
      --temp $TEMP --nvt-ns 0.2 --npt-ns 1.0 2>&1 | tail -3
  python "$HERE/gen_gamd_config.py" --prmtop "$PRM" \
      --rst7 "$RUNS/$TAG/sys1_${drug}_equil.rst7" \
      --outdir "$RUNS/$TAG/gamd_out" --out "$RUNS/$TAG/config.xml" \
      --prod-ns $PROD_S1 --no-min --boost-type lower-dual --seed $SEED \
      --report-ps 50 --temp $TEMP 2>&1 | tail -2
  echo "  [$TAG] built + equilibrated + config ready"
done

# launch S1 on GPUs 0/1/2
i=0
for drug in niraparib olaparib rucaparib; do
  TAG="sys1_${drug}"
  rm -f "$RUNS/$TAG/gamd_out/gamd.log" "$RUNS/$TAG/gamd_out/output.dcd"
  CUDA_VISIBLE_DEVICES=$i nohup "$ROOT/tools/gamd-openmm/gamdRunner" xml \
      "$RUNS/$TAG/config.xml" > "$RUNS/$TAG/runner.log" 2>&1 < /dev/null &
  echo "  [$TAG] gamdRunner PID $! -> GPU $i"
  i=$((i+1))
  sleep 5
done
while pgrep -f gamdRunner > /dev/null; do sleep 600; done
for drug in niraparib olaparib rucaparib; do
  TAG="sys1_${drug}"
  cp "$RUNS/$TAG/gamd_out/gamd.log" "$ARCHIVE/sys1_${drug}_gamd.log"
  echo "  [$TAG] gamd.log archived ($(wc -l < "$ARCHIVE/sys1_${drug}_gamd.log") lines)"
done

echo "==== ALL 6 SYSTEMS DONE. Logs in $ARCHIVE ===="
ls -la "$ARCHIVE"
