#!/usr/bin/env bash
# deploy_s1_only.sh — run only the S1 batch (niraparib/olaparib/rucaparib,
# 24 ns, full rebuild). The S2 batch already completed and its logs are
# archived; this script must NOT touch S2.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REC_S1="$ROOT/docking/receptors/6VKK_A_receptor.pdb"
LIG_S1="$ROOT/docking/results"
RUNS="$ROOT/runs_rerun"
ARCHIVE="$ROOT/results/rerun_logs"
mkdir -p "$RUNS" "$ARCHIVE"

PROD_S1=24
SEED=42
TEMP=300

echo "==== [S1] niraparib/olaparib/rucaparib (24 ns, full rebuild) ===="
for drug in niraparib olaparib rucaparib; do
  TAG="sys1_${drug}"
  PRM="$RUNS/${TAG}.prmtop"
  INC="$RUNS/${TAG}.inpcrd"
  if [ ! -f "$PRM" ]; then
    python "$HERE/build_system.py" --receptor "$REC_S1" \
        --ligand "$LIG_S1/${drug}_docked.fixed.sdf" --out "$RUNS/$TAG" 2>&1 | tail -4
  fi
  EQ="$RUNS/$TAG/${TAG}_equil.rst7"
  if [ ! -f "$EQ" ]; then
    python "$HERE/equilibrate.py" --prmtop "$PRM" --inpcrd "$INC" \
        --out "$RUNS/$TAG/${TAG}_equil" 2>&1 | tail -3
  fi
  python "$HERE/gen_gamd_config.py" --prmtop "$PRM" --rst7 "$EQ" \
      --outdir "$RUNS/$TAG/gamd_out" --out "$RUNS/$TAG/config.xml" \
      --prod-ns $PROD_S1 --no-min --boost-type lower-dual --seed $SEED \
      --report-ps 50 --temp $TEMP 2>&1 | tail -2
  echo "  [$TAG] config ready"
done

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
echo "==== S1 batch launched; waiting for completion ===="
while pgrep -f gamdRunner > /dev/null; do sleep 300; done
for drug in niraparib olaparib rucaparib; do
  TAG="sys1_${drug}"
  cp "$RUNS/$TAG/gamd_out/gamd.log" "$ARCHIVE/sys1_${drug}_gamd.log"
  echo "  [$TAG] gamd.log archived ($(wc -l < "$ARCHIVE/sys1_${drug}_gamd.log") lines)"
done
echo "S1_DONE"
