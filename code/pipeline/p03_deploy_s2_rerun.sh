#!/usr/bin/env bash
# p03_deploy_s2_rerun.sh — catch-up S2 runs for AZD5305 + veliparib
#
# Purpose:   Re-run the AZD5305 + veliparib S2 systems (the legacy runs had
#            only 2.6 ns) to 26 ns production, matching the other legacy S2
#            systems. Run on the same workstation after p02 finishes.
# GPU:       AZD5305=GPU0, veliparib=GPU1 (300 K, consistent with legacy S2).
# Author:    Tao Zhu (tjogzt@gmail.com)
# Created:   2026-09-17
# Inputs:    ../../data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb
#            ../../data/01_modeling/docking_system2/mol2/<drug>_best.mol2
# Outputs:   ../../runs_s2/sys2_<drug>/{build, equil, config.xml, prod logs}
# Run:       bash code/pipeline/p03_deploy_s2_rerun.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REC="$HERE/../../data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb"
LIGDIR="$HERE/../../data/01_modeling/docking_system2"
RUNS="$HERE/../../runs_s2"
mkdir -p "$RUNS"
PROD_NS=26
TEMP=300

declare -A GPU=( [AZD5305]=0 [veliparib]=1 )

for drug in AZD5305 veliparib; do
  TAG="sys2_${drug}"
  PRM="$RUNS/$TAG/build/${drug}_best.prmtop"
  INC="$RUNS/$TAG/build/${drug}_best.inpcrd"

  if [ -f "$PRM" ] && [ -f "$INC" ]; then
    echo "==== [$TAG] 1/4 build skipped (prmtop/inpcrd already present) ===="
  else
    echo "==== [$TAG] 1/4 tleap build ===="
    python "$HERE/p01_build_system2_tleap.py" --receptor "$REC" \
        --ligand "$LIGDIR/mol2/${drug}_best.mol2" --out "$RUNS/$TAG" 2>&1 | tail -6
  fi

  echo "==== [$TAG] 2/4 equilibration (NVT 0.2ns + NPT 1.0ns, ${TEMP}K) ===="
  CUDA_VISIBLE_DEVICES=${GPU[$drug]} python "$HERE/p04_equilibrate.py" \
      --prmtop "$PRM" --inpcrd "$INC" \
      --out "$RUNS/$TAG/${TAG}_equil" --temp $TEMP --nvt-ns 0.2 --npt-ns 1.0 2>&1 | tail -4

  echo "==== [$TAG] 3/4 GaMD config (${PROD_NS}ns production, 50ps frames, lower-dual, ${TEMP}K) ===="
  python "$HERE/p05_gen_gamd_config.py" --prmtop "$PRM" \
      --rst7 "$RUNS/$TAG/${TAG}_equil.rst7" --outdir "$RUNS/$TAG/gamd_out" \
      --out "$RUNS/$TAG/config.xml" --prod-ns $PROD_NS --no-min \
      --boost-type lower-dual --seed 42 --report-ps 50 --temp $TEMP 2>&1 | tail -3

  echo "==== [$TAG] 4/4 launching GaMD (GPU ${GPU[$drug]}) ===="
  rm -rf "$RUNS/$TAG/gamd_out/gamd.log" "$RUNS/$TAG/gamd_out/output.dcd" \
         "$RUNS/$TAG/gamd_out/gamd_restart.checkpoint" 2>/dev/null || true
  CUDA_VISIBLE_DEVICES=${GPU[$drug]} nohup "$HERE/../../tools/gamd-openmm/gamdRunner" xml \
      "$RUNS/$TAG/config.xml" > "$RUNS/$TAG/gamd.log" 2>&1 < /dev/null &
  echo "  [$TAG] gamdRunner PID $! → GPU ${GPU[$drug]}"
  sleep 5
done

echo "==== both catch-up runs launched ===="
sleep 15
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader
