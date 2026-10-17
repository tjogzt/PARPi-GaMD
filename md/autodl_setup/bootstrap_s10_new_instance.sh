#!/bin/bash
# bootstrap_s10_new_instance.sh — S10 前置就绪（在"换网点新建"的实例上执行）
#
# 适用：无法跨网点克隆，新开干净实例后从零恢复 S10 运行环境。
# 前提：① 已把 s10_bundle.tgz 上传到 /root/s10_bundle.tgz（内含两个 niraparib S2
#         体系的 prmtop/restart/config/gamd.log）；
#       ② 实例有外网（git clone + conda 环境创建需要）。
# 用法：bash bootstrap_s10_new_instance.sh   （任意 cwd）
# 之后：bash /root/autodl-tmp/PARPi_design/md/autodl_setup/deploy_s2_nira_multiseed_v2.sh
set -euo pipefail

R=/root/autodl-tmp/PARPi_design
BUNDLE=/root/s10_bundle.tgz
TAG=v2.2
REPO=https://github.com/tjogzt/PARPi-GaMD.git

echo "=== [1/6] 代码 ==="
if [ -d "$R" ]; then
  echo "exists: $R (skip clone)"
else
  git clone --depth 1 --branch "$TAG" "$REPO" "$R"
fi
cd "$R"

echo "=== [2/6] conda 环境 ==="
CONDA_BASE="${CONDA_BASE:-/root/miniconda3}"
if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  echo "ERROR: conda not found at $CONDA_BASE（AutoDL 镜像一般自带 /root/miniconda3）"
  exit 1
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
if conda env list | awk '{print $1}' | grep -qx gamd; then
  echo "env gamd exists"
else
  conda env create -f "$R/md/autodl_setup/environment_gamd_min.yml"
fi
conda activate gamd
python -V
python -c "import openmm, parmed, numpy, scipy, mdtraj; print('openmm', openmm.__version__, '| parmed', parmed.__version__)"

echo "=== [3/6] 解包 S10 体系文件 ==="
[ -f "$BUNDLE" ] || { echo "ERROR: $BUNDLE 不存在，请先上传 s10_bundle.tgz"; exit 1; }
tar xzf "$BUNDLE" -C "$R"
echo "--- 关键文件 ---"
ls -lh "$R/runs/sys2_niraparib/sys2_niraparib.prmtop" \
       "$R/runs/sys2_niraparib/sys2_niraparib_start.rst7" \
       "$R/runs_nira_prot/sys2_niraparib_prot/sys2_niraparib_prot_equil.rst7" \
       "$R/runs_nira_prot/sys2_niraparib_prot/build/niraparib_best_prot.prmtop"

echo "=== [4/6] gamdRunner 冒烟测试（-h，无需 GPU）==="
"$R/tools/gamd-openmm/gamdRunner" -h > /dev/null 2>&1 && echo "gamdRunner OK" || { echo "gamdRunner FAILED"; exit 1; }

echo "=== [5/6] GPU 与磁盘 ==="
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader || echo "WARN: no GPU visible"
df -h /root/autodl-tmp | tail -1
NGC=$(nvidia-smi -L 2>/dev/null | wc -l || echo 0)
echo "GPU count = $NGC"
[ "$NGC" -ge 1 ] || { echo "ERROR: 无 GPU，开正式模式后再跑"; exit 1; }

echo "=== [6/6] 就绪 ==="
echo "下一步：bash $R/md/autodl_setup/deploy_s2_nira_multiseed_v2.sh"
