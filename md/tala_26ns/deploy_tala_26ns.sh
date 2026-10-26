#!/bin/bash
# deploy_tala_26ns.sh — launch a 26-ns talazoparib S2 production run
# from the archived 12-ns equilibrated restart, same GaMD protocol as the panel.
# Idempotent: if gamd_out/output.dcd already exists, it refuses to overwrite.
set -euo pipefail

ROOT=/root/autodl-tmp/PARPi_design/runs_extension/tala_s2_26ns
cd "$ROOT"

CONDA_BASE="${CONDA_BASE:-}"
if [ -z "$CONDA_BASE" ]; then
  for c in /root/miniconda3 /root/autodl-tmp/miniconda3 /root/autodl-tmp/conda; do
    if [ -f "$c/etc/profile.d/conda.sh" ]; then CONDA_BASE="$c"; break; fi
  done
fi
[ -n "$CONDA_BASE" ] || { echo "FATAL: conda not found"; exit 1; }
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate gamd
PY="$CONDA_PREFIX/bin/python"

echo "=== [0] environment ==="
$PY -c 'import openmm; print("openmm", openmm.__version__); print("platforms", [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])'
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
NGC=$(nvidia-smi -L | wc -l)
[ "$NGC" -ge 1 ] || { echo "FATAL: no GPU"; exit 1; }

echo "=== [1] inputs ==="
for f in tala.prmtop tala_s2_equil12ns.rst7 config_tala_s2_26ns_seed2026.xml; do
  [ -s "$f" ] || { echo "FATAL missing $f"; exit 1; }
done
$PY - <<'PY'
import parmed as pmd
t = pmd.load_file('tala.prmtop'); r = pmd.load_file('tala_s2_equil12ns.rst7')
print('atoms prmtop/rst7:', len(t.atoms), r.natom)
assert len(t.atoms) == r.natom == 285796, 'atom-count mismatch'
print('net charge:', sum(a.charge for a in t.atoms))
PY

mkdir -p gamd_out
RESUME_FLAG=""
if [ -s gamd_out/output.dcd ]; then
  if [ "${RESUME:-0}" = "1" ] && [ -s gamd_out/gamd_restart.checkpoint ]; then
    echo "RESUME=1: appending from checkpoint."
    RESUME_FLAG="-r"
  else
    echo "EXISTING gamd_out/output.dcd found; refusing to overwrite (set RESUME=1 to continue)."
    ls -lh gamd_out/output.dcd
    exit 0
  fi
fi

echo "=== [2] launch seed 2026, 26 ns production ${RESUME_FLAG} ==="
export CUDA_VISIBLE_DEVICES=0
nohup "$PY" tools/gamdRunner xml config_tala_s2_26ns_seed2026.xml $RESUME_FLAG \
  > run_seed2026.log 2>&1 </dev/null &
echo "PID=$!"
sleep 5
ps -o pid,etime,pcpu,cmd -p "$!" || true
echo "=== launched; monitor gamd_out/gamd.log and output.dcd ==="
