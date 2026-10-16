#!/bin/bash
# run_review_replicates.sh — REV-1/REV-2 compute batch (run on the AutoDL instance)
# Targets:
#   REV-2: talazoparib S1 main-pipeline replicates 3x200 ns (seed 49/7/123)
#   REV-1: veliparib   S1 extension + replicates 2x400 ns (seed 49/7)
#   Analysis: main-pipeline CV (HD 0-125 / ART 126-344, internal numbering) + PyReweighting C3
# Usage: bash run_review_replicates.sh  (cwd = /root/autodl-tmp/PARPi_design)
set -euo pipefail
R=/root/autodl-tmp/PARPi_design
cd "$R"
PY=/usr/bin/python3
export CUDA_VISIBLE_DEVICES=""

echo "=== [0] env check ==="
$PY -c 'import openmm; print([p.getName() for p in [openmm.Platform.getPlatform(i) for i in range(openmm.Platform.getNumPlatforms())]])'
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
df -h /root/autodl-tmp | tail -1

build_one() {  # $1=TAG $2=ligand_sdf (empty = apo) $3=prod_ns
  local TAG=$1 LIG=${2:-} PRODNS=$3
  local SYS=$R/runs/$TAG
  mkdir -p "$SYS"
  # --- build ---
  if [ ! -f "$SYS/${TAG}.prmtop" ]; then
    echo "=== build $TAG ==="
    if [ -n "$LIG" ]; then
      $PY md/autodl_setup/build_system.py --receptor docking/receptors/6VKK_A_receptor.pdb \
        --ligand "$LIG" --out "$SYS" 2>&1 | tail -3
    fi
  fi
  # --- equilibrate (cMD 2ns + GaMD equil 10ns) ---
  if [ ! -f "$SYS/${TAG}_equil.rst7" ]; then
    echo "=== equilibrate $TAG ==="
    $PY md/autodl_setup/equilibrate.py --prmtop "$SYS/${TAG}.prmtop" --out "$SYS" 2>&1 | tail -3
  fi
}

launch_prod() {  # $1=TAG $2=seed $3=prod_ns $4=gpu
  local TAG=$1 SEED=$2 PRODNS=$3 GPU=$4
  local SYS=$R/runs/$TAG
  local CFG=$SYS/config_${TAG}_s${SEED}.xml
  $PY md/autodl_setup/gen_gamd_config.py --prmtop "$SYS/${TAG}.prmtop" \
    --rst7 "$SYS/${TAG}_equil.rst7" --outdir "$SYS/gamd_out_s${SEED}" \
    --out "$CFG" --prod-ns "$PRODNS" --no-min --boost-type lower-dual \
    --seed "$SEED" --report-ps 50
  sed -i "s|_equil.rst7|${TAG}_equil.rst7|g" "$CFG"
  grep -q "${TAG}_equil.rst7" "$CFG" || { echo "RST7 PATH BUG in $CFG"; exit 1; }
  echo "=== launch $TAG seed=$SEED ${PRODNS}ns GPU=$GPU ==="
  CUDA_VISIBLE_DEVICES=$GPU nohup "$R/tools/gamd-openmm/gamdRunner" xml "$CFG" \
    > /tmp/gamd_${TAG}_s${SEED}.log 2>&1 &
  sleep 20
  ps aux | grep -c "[g]amdRunner" || true
}

echo "=== [1] build talazoparib S1 ==="
build_one sys1_talazoparib docking/results/talazoparib_docked.fixed.sdf 200
echo "=== [2] build veliparib S1 ==="
build_one sys1_veliparib docking/results/veliparib_docked.sdf 400

echo "=== [3] launch productions (2 jobs/GPU max) ==="
# GPU0: tala s49 (200) -> tala s7 (200); GPU1: veli s49 (400) -> veli s7 (400); tala s123 queued
launch_prod sys1_talazoparib 49 200 0
launch_prod sys1_veliparib 49 400 1

echo "=== [4] queue remaining (launch after earlier ones finish; run queue_monitor.sh) ==="
cat > /tmp/queue_review.sh <<'QEOF'
#!/bin/bash
# launch the follow-up seeds once s49 completes
R=/root/autodl-tmp/PARPi_design
while true; do
  for pair in "sys1_talazoparib:7:200:0:49" "sys1_veliparib:7:400:1:49" "sys1_talazoparib:123:200:0:7"; do
    IFS=: read TAG SEED NS GPU DONE_S <<< "$pair"
    SYS=$R/runs/$TAG
    if [ -f "$SYS/gamd_out_s${SEED}/gamd.log" ]; then continue; fi
    DONEFILE=$SYS/gamd_out_s${DONE_S}/gamd.log
    if [ -f "$DONEFILE" ]; then
      STEP=$(tail -1 "$DONEFILE" | awk '{print $2}')
      if [ "$STEP" -ge 6000000 ]; then
        CFG=$SYS/config_${TAG}_s${SEED}.xml
        if [ ! -f "$CFG" ]; then
          /usr/bin/python3 md/autodl_setup/gen_gamd_config.py --prmtop "$SYS/${TAG}.prmtop" \
            --rst7 "$SYS/${TAG}_equil.rst7" --outdir "$SYS/gamd_out_s${SEED}" \
            --out "$CFG" --prod-ns "$NS" --no-min --boost-type lower-dual \
            --seed "$SEED" --report-ps 50
          sed -i "s|_equil.rst7|${TAG}_equil.rst7|g" "$CFG"
        fi
        echo "launch $TAG seed=$SEED"
        CUDA_VISIBLE_DEVICES=$GPU nohup $R/tools/gamd-openmm/gamdRunner xml "$CFG" > /tmp/gamd_${TAG}_s${SEED}.log 2>&1 &
        sleep 20
      fi
    fi
  done
  sleep 600
done
QEOF
nohup bash /tmp/queue_review.sh > /tmp/queue_review.log 2>&1 &
echo "DEPLOY DONE. Productions launching; queue monitor running."
