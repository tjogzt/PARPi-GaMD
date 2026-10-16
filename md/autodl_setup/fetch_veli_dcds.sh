#!/bin/bash
# fetch_veli_dcds.sh — fetch veli_prot trajectories (sys2 2.6G + sys1 0.5G) to the data volume.
# Uses rsync (resumable) when available, else scp with retries. Runs from anywhere.
set -u
DEST=/Volumes/tjogzt4T/PARPi_data/veli_prot
H="root@connect.nmb2.seetacloud.com"
P=38683
SSHOPTS="-p $P -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o ConnectTimeout=25"
SCPOPTS="-P $P -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o ConnectTimeout=25"
mkdir -p "$DEST/sys2_veliparib_prot" "$DEST/sys1_veliparib_prot"

fetch() {  # fetch <remote_path> <local_dir> <name>
  local rp="$1" ld="$2" nm="$3"
  if command -v rsync >/dev/null && ssh $SSHOPTS "$H" "command -v rsync" >/dev/null 2>&1; then
    rsync -avP --partial -e "ssh $SSHOPTS" "$H:$rp" "$ld/$nm" && return 0
  fi
  for i in 1 2 3; do
    scp $SCPOPTS "$H:$rp" "$ld/$nm" && return 0
    echo "retry $i for $nm"; sleep 30
  done
  return 1
}

echo "== sys2 dcd =="
fetch "/root/autodl-tmp/PARPi_design/runs_veli_prot/sys2_veliparib_prot/gamd_out/output.dcd" "$DEST/sys2_veliparib_prot" "output.dcd"
echo "== sys1 dcd =="
fetch "/root/autodl-tmp/PARPi_design/runs_veli_prot/sys1_veliparib_prot/gamd_out/output.dcd" "$DEST/sys1_veliparib_prot" "output.dcd"
echo "== exts =="
for d in sys2_veliparib_prot sys1_veliparib_prot; do
  ls -la "$DEST/$d/output.dcd" 2>/dev/null | awk '{print $5, $9}'
done
echo "FETCH_DONE"
