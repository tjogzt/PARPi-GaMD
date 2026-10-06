#!/usr/bin/env bash
# p02_deploy_s2.sh — deploy the three extension inhibitors' S2 GaMD runs
#
# Purpose:   One-shot build -> equilibrate -> GaMD launch for the three
#            extension-inhibitor S2 systems (full-length PARP1-DNA) on a 3-GPU
#            workstation; each drug gets its own GPU (fluzoparib=GPU0,
#            pamiparib=GPU1, senaparib=GPU2).
# Protocol:  300 K (consistent with the legacy S1/S2 runs; legacy S2 weights
#            beta=1.677571 -> 300.15 K verified), 22 ns production, 50 ps
#            frames, lower-dual boost, seed 42.
# Created:   2026-09-15
# Inputs:    ../../data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb
#            ../../data/01_modeling/docking_system2/mol2/<drug>_best.mol2
#            sibling scripts: p01 build, p04 equilibrate, p05 gen_gamd_config
# Outputs:   ../../runs_s2/sys2_<drug>/{build, equil, config.xml, prod logs}
# Run:       bash code/pipeline/p02_deploy_s2.sh   (GPU workstation, env 'gamd')
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REC="$HERE/../../data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb"
LIGDIR="$HERE/../../data/01_modeling/docking_system2"
RUNS="$HERE/../../runs_s2"
mkdir -p "$RUNS"
PROD_NS=22
TEMP=300

declare -A GPU=( [fluzoparib]=0 [pamiparib]=1 [senaparib]=2 )

for drug in fluzoparib pamiparib senaparib; do
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

echo "==== all three S2 runs launched ===="
sleep 15
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader
