#!/bin/bash
# ============================================================================
# GaMD Replicate Trajectories for PARP1 Trapping Paper
# Run on AutoDL: bash replicate_runs.sh
# 
# Priority systems (S1 CAT-only, ~30K atoms, ~100 ns each):
#   1. talazoparib (Type II extreme, strongest trapper)   — 2 repeat
#   2. veliparib   (Type III extreme, weakest trapper)     — 2 repeat + extend to 300 ns
#   3. AZD5305     (PARP1-selective, blind prediction)     — 2 repeat
#
# Expected GPU time: ~2-3 days on 2×RTX 4090
# ============================================================================
set -e
export LANG=C

CONDA_ENV="gamd"
WORKDIR="/root/autodl-tmp/PARPi_design/runs"
SCRIPT_DIR="/root/autodl-tmp/PARPi_design/md/autodl_setup"

echo "========================================"
echo " PARP1 GaMD REPLICATE TRAJECTORIES"
echo "========================================"
echo "Instance: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo ""

source activate $CONDA_ENV 2>/dev/null || conda activate $CONDA_ENV

# Define replicate runs: (system, ligand, seed, ns)
declare -A RUNS
RUNS=(
  "sys1_talazoparib_rep1|sys1|talazoparib|42|100"
  "sys1_talazoparib_rep2|sys1|talazoparib|123|100"
  "sys1_veliparib_rep1|sys1|veliparib|42|100"
  "sys1_veliparib_rep2|sys1|veliparib|123|100"
  "sys1_AZD5305_rep1|sys1|AZD5305|42|100"
  "sys1_AZD5305_rep2|sys1|AZD5305|123|100"
)

for run_info in "${RUNS[@]}"; do
  IFS='|' read -r NAME SYSTEM LIGAND SEED NS <<< "$run_info"
  
  RUNDIR="$WORKDIR/${NAME}"
  if [ -d "$RUNDIR" ] && [ -f "$RUNDIR/done.tag" ]; then
    echo "SKIP $NAME (already complete)"
    continue
  fi
  
  echo "=== Submitting: $NAME (seed=$SEED, ${NS}ns) ==="
  mkdir -p "$RUNDIR"
  
  # Build system with new seed
  python3 $SCRIPT_DIR/build_system.py \
    --system $SYSTEM --ligand $LIGAND \
    --seed $SEED --output $RUNDIR
  
  # Equilibrate
  python3 $SCRIPT_DIR/equilibrate.py \
    --input $RUNDIR/system.prmtop \
    --coords $RUNDIR/system.inpcrd \
    --seed $SEED
  
  # Run GaMD production
  python3 $SCRIPT_DIR/run_prod.sh \
    --input $RUNDIR/system.prmtop \
    --coords $RUNDIR/equil.rst7 \
    --ns $NS --seed $SEED \
    --output $RUNDIR
  
  # Post-run analysis (PMF reweighting)
  python3 $SCRIPT_DIR/analyze_gamd.py \
    --trajectory $RUNDIR/prod.nc \
    --topology $RUNDIR/system.prmtop \
    --cv HD_ART_distance \
    --output $RUNDIR/analysis
  
  touch "$RUNDIR/done.tag"
  echo "COMPLETE $NAME"
  echo ""
done

echo ""
echo "========================================"
echo " ALL REPLICATES COMPLETE"
echo "========================================"
echo ""
echo "Download analysis results:"
echo "  scp -r root@<autodl-ip>:/root/autodl-tmp/PARPi_design/runs/sys1_*_rep*/analysis ./results/replicates/"
