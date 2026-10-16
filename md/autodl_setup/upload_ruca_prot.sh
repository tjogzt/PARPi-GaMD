#!/usr/bin/env bash
# upload_ruca_prot.sh — upload rucaparib protonated rerun inputs to AutoDL (nmb1)
# Run locally. Target instance: ssh <alias> (host/port kept out of the repository; configure in ~/.ssh/config).
set -euo pipefail
LOCAL_ROOT="${LOCAL_ROOT:-$HOME/clacky_workspace/PARPi_design}"
REMOTE="autodl-nmb1"
RDIR="/root/autodl-tmp/PARPi_design"

scp -q "$LOCAL_ROOT/docking/results_system2/mol2/rucaparib_best_prot.mol2" \
  "$REMOTE:$RDIR/docking/results_system2/mol2/"
scp -q "$LOCAL_ROOT/docking/results/rucaparib_docked_prot.sdf" \
  "$REMOTE:$RDIR/docking/results/"
scp -q "$LOCAL_ROOT/md/autodl_setup/deploy_ruca_prot.sh" \
  "$REMOTE:$RDIR/md/autodl_setup/"
scp -q "$LOCAL_ROOT/scripts/protonate_rucaparib.py" \
  "$REMOTE:$RDIR/scripts/"
echo "upload done; next on instance:"
echo "  ssh $REMOTE 'bash /root/autodl-tmp/PARPi_design/md/autodl_setup/deploy_ruca_prot.sh'"
