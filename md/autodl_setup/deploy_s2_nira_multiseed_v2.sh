#!/bin/bash
# deploy_s2_nira_multiseed_v2.sh — S10 补充计算批（换网点新建实例版）
# 与 v1 的差异：① conda 环境路径可配置（默认 /root/miniconda3，AutoDL 镜像自带）；
#   ② 体系文件发现兼容三种布局（*equil*.rst7 或 *_start.rst7；顶层或 build/ 下的 prmtop）；
#   ③ GPU 列表按 nvidia-smi 实际卡数自动生成（默认 1 job/卡；可显式覆盖）。
# 背景：niraparib 双态 S2 CV1 多种子复制（判定单轨迹值的跨种子稳定性）。
# 目标：sys2_niraparib（中性）与 sys2_niraparib_prot（+1）各补 2 个独立 seed。
# 用法：bash deploy_s2_nira_multiseed_v2.sh   （任意 cwd）
set -euo pipefail
R=/root/autodl-tmp/PARPi_design
cd "$R"
CONDA_BASE="${CONDA_BASE:-/root/miniconda3}"
source "$CONDA_BASE/etc/profile.d/conda.sh" && conda activate gamd
PY="$CONDA_PREFIX/bin/python"

PRODNS=${PRODNS:-26}           # 26 ns，与面板标准生产长度一致
SEEDS_NEW="${SEEDS_NEW:-7 123}" # 新增 seed（既有中性/质子化重跑均用 seed 42，不冲突）
declare -A SYS=( [niraparib]=runs/sys2_niraparib [niraparib_prot]=runs_nira_prot/sys2_niraparib_prot )

echo "=== [0] 环境与既有 seed 核验 ==="
$PY -c 'import openmm; print("platforms:", [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])'
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
NGC=$(nvidia-smi -L | wc -l)
if [ -n "${GPU_LIST:-}" ]; then
  read -r -a GPUS <<< "$GPU_LIST"
else
  GPUS=()
  for ((g=0; g<NGC; g++)); do GPUS+=("$g"); done   # 1 job/卡（最稳；如需 2/卡自行改）
fi
echo "GPU_LIST = ${GPUS[*]}  (per-card jobs = 1)"
for st in niraparib niraparib_prot; do
  d="$R/${SYS[$st]}"
  echo "-- $st @ $d"
  [ -d "$d" ] || { echo "MISSING $d"; exit 1; }
  grep -o "<random-seed>[0-9]*</random-seed>" "$d"/config.xml 2>/dev/null | head -1 || echo "  (no config seed found)"
done

launch_one() {  # $1=state $2=seed $3=gpu
  local ST=$1 SEED=$2 GPU=$3
  local SYS_DIR="$R/${SYS[$ST]}"
  local TAG="sys2_niraparib${ST/niraparib/}"
  local OUT="$SYS_DIR/gamd_out_s${SEED}_r14"
  local CFG="$SYS_DIR/config_${TAG}_s${SEED}_r14.xml"
  [ -d "$OUT" ] && { echo "SKIP existing $OUT"; return 0; }
  # 起点构象：优先 equil，其次 *_start（两种历史命名）
  local EQUIL PRMTOP
  EQUIL=$(ls "$SYS_DIR"/*equil*.rst7 2>/dev/null | head -1 || true)
  [ -n "${EQUIL:-}" ] || EQUIL=$(ls "$SYS_DIR"/*_start.rst7 2>/dev/null | head -1 || true)
  PRMTOP=$(ls "$SYS_DIR"/*.prmtop 2>/dev/null | head -1 || true)
  [ -n "${PRMTOP:-}" ] || PRMTOP=$(ls "$SYS_DIR"/build/*.prmtop 2>/dev/null | head -1 || true)
  [ -n "${EQUIL:-}" ] && [ -n "${PRMTOP:-}" ] || { echo "MISSING prmtop/restart in $SYS_DIR"; exit 1; }
  echo "  $TAG: prmtop=$(basename "$PRMTOP")  start=$(basename "$EQUIL")"
  $PY md/autodl_setup/gen_gamd_config.py --prmtop "$PRMTOP" --rst7 "$EQUIL" \
      --outdir "$OUT" --out "$CFG" --prod-ns "$PRODNS" --no-min --boost-type lower-dual \
      --seed "$SEED" --report-ps 50
  # 绝对路径整元素重写（防双重拼接）
  $PY - "$CFG" "$PRMTOP" "$EQUIL" <<'EOP'
import re, sys
cfg, prm, rst = sys.argv[1:4]
x = open(cfg).read()
x = re.sub(r'<coordinates type="rst7">.*?</coordinates>',
           f'<coordinates type="rst7">{rst}</coordinates>', x)
x = re.sub(r'<topology>.*?</topology>', f'<topology>{prm}</topology>', x)
open(cfg, "w").write(x)
EOP
  echo "=== launch $TAG seed=$SEED ${PRODNS}ns GPU=$GPU ==="
  CUDA_VISIBLE_DEVICES=$GPU nohup "$R/tools/gamd-openmm/gamdRunner" xml "$CFG" \
      > "/tmp/gamd_${TAG}_s${SEED}_r14.log" 2>&1 </dev/null &
  sleep 2
}

i=0
for SEED in $SEEDS_NEW; do
  for ST in niraparib niraparib_prot; do
    launch_one "$ST" "$SEED" "${GPUS[$((i % ${#GPUS[@]}))]}"
    i=$((i+1))
  done
done

echo "=== 队列完成（后台运行）==="
ps -o pid,etime,cmd -C python3 | grep -i gamd | head -6 || true
echo "监控：nvidia-smi + DCD 帧数增速 + gamd_out_s*/gamd.log 步数；完成判据 = 生产步数达标（6.0e6 + PRODNS×5e5）"
