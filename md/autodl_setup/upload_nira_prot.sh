#!/usr/bin/env bash
# upload_nira_prot.sh — upload niraparib protonated rerun inputs to AutoDL
# Run locally. Target instance: ssh <alias> (host/port kept out of the repository; configure in ~/.ssh/config).
set -euo pipefail
LOCAL_ROOT="${LOCAL_ROOT:-$HOME/clacky_workspace/PARPi_design}"
REMOTE="autodl"
RDIR="/root/autodl-tmp/PARPi_design"

scp -q \
  "$LOCAL_ROOT/docking/results_system2/mol2/niraparib_best_prot.mol2" \
  "$REMOTE:$RDIR/docking/results_system2/mol2/"
scp -q \
  "$LOCAL_ROOT/docking/results/niraparib_docked_prot.sdf" \
  "$REMOTE:$RDIR/docking/results/"
scp -q \
  "$LOCAL_ROOT/data/01_modeling/parp1_dna_complex/PARPi_full_DNA_Zn.pdb" \
  "$REMOTE:$RDIR/data/01_modeling/parp1_dna_complex/"
scp -q \
  "$LOCAL_ROOT/docking/receptors/6VKK_A_receptor.pdb" \
  "$REMOTE:$RDIR/docking/receptors/"
scp -q \
  "$LOCAL_ROOT/md/autodl_setup/deploy_nira_prot.sh" \
  "$REMOTE:$RDIR/md/autodl_setup/"
echo "upload done; next on instance:"
echo "  ssh $REMOTE 'bash /root/autodl-tmp/PARPi_design/md/autodl_setup/deploy_nira_prot.sh'"
