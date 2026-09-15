#!/usr/bin/env bash
# S2 部署: 三个新药的全长 PARP1-DNA 体系 build→equil→GaMD 一条龙
# 用法: bash deploy_s2.sh   (在 AutoDL 3×4090 实例, gamd 环境已激活)
# 每药独立 GPU: fluzoparib=GPU0, pamiparib=GPU1, senaparib=GPU2
# 协议: 300K (与旧 S1/S2 一致, 旧 S2 权重 β=1.677571→300.15K 已核实)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REC="$HERE/../data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb"
LIGDIR="$HERE/../docking/results_system2"
RUNS="$HERE/../runs_s2"
mkdir -p "$RUNS"
PROD_NS=22
TEMP=300

declare -A GPU=( [fluzoparib]=0 [pamiparib]=1 [senaparib]=2 )

for drug in fluzoparib pamiparib senaparib; do
  TAG="sys2_${drug}"
  PRM="$RUNS/$TAG/build/${drug}_best.prmtop"
  INC="$RUNS/$TAG/build/${drug}_best.inpcrd"

  if [ -f "$PRM" ] && [ -f "$INC" ]; then
    echo "==== [$TAG] 1/4 构建跳过 (prmtop/inpcrd 已存在) ===="
  else
    echo "==== [$TAG] 1/4 tleap 构建 ===="
    python "$HERE/build_system2_tleap.py" --receptor "$REC" \
        --ligand "$LIGDIR/mol2/${drug}_best.mol2" --out "$RUNS/$TAG" 2>&1 | tail -6
  fi

  echo "==== [$TAG] 2/4 平衡 (NVT 0.2ns + NPT 1.0ns, ${TEMP}K) ===="
  CUDA_VISIBLE_DEVICES=${GPU[$drug]} python "$HERE/equilibrate.py" \
      --prmtop "$PRM" --inpcrd "$INC" \
      --out "$RUNS/$TAG/${TAG}_equil" --temp $TEMP --nvt-ns 0.2 --npt-ns 1.0 2>&1 | tail -4

  echo "==== [$TAG] 3/4 GaMD config (${PROD_NS}ns 生产, 50ps 帧率, lower-dual, ${TEMP}K) ===="
  python "$HERE/gen_gamd_config.py" --prmtop "$PRM" \
      --rst7 "$RUNS/$TAG/${TAG}_equil.rst7" --outdir "$RUNS/$TAG/gamd_out" \
      --out "$RUNS/$TAG/config.xml" --prod-ns $PROD_NS --no-min \
      --boost-type lower-dual --seed 42 --report-ps 50 --temp $TEMP 2>&1 | tail -3

  echo "==== [$TAG] 4/4 启动 GaMD (GPU ${GPU[$drug]}) ===="
  rm -rf "$RUNS/$TAG/gamd_out/gamd.log" "$RUNS/$TAG/gamd_out/output.dcd" \
         "$RUNS/$TAG/gamd_out/gamd_restart.checkpoint" 2>/dev/null || true
  CUDA_VISIBLE_DEVICES=${GPU[$drug]} nohup "$HERE/../tools/gamd-openmm/gamdRunner" xml \
      "$RUNS/$TAG/config.xml" > "$RUNS/$TAG/gamd.log" 2>&1 < /dev/null &
  echo "  [$TAG] gamdRunner PID $! → GPU ${GPU[$drug]}"
  sleep 5
done

echo "==== 三个 S2 全部启动 ===="
sleep 15
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader
