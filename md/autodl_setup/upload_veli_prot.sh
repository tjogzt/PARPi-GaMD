#!/usr/bin/env bash
# upload_veli_prot.sh — upload the inputs required by the veliparib protonated rerun
# Usage: bash md/autodl_setup/upload_veli_prot.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
LOCAL_ROOT="$(cd "$HERE/../.." && pwd)"
REMOTE="root@connect.nmb1.seetacloud.com"
PORT=30694
RDIR="/root/autodl-tmp/PARPi_design"

scp -P "$PORT" "$LOCAL_ROOT/docking/results_system2/mol2/veliparib_best_prot.mol2" \
    "$REMOTE:$RDIR/docking/results_system2/mol2/"
scp -P "$PORT" "$LOCAL_ROOT/docking/results/veliparib_docked_prot.sdf" \
    "$REMOTE:$RDIR/docking/results/"
scp -P "$PORT" "$LOCAL_ROOT/md/autodl_setup/deploy_veli_prot.sh" \
    "$REMOTE:$RDIR/md/autodl_setup/"
echo "UPLOAD OK"
