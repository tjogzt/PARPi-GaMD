#!/bin/bash
# pull_ruca_prot.sh — download rucaparib protonated rerun products (sys1|sys2|both)
# Run locally (Mac): bash pull_ruca_prot.sh sys1
# Layout mirrors runs_nira_prot: <dest>/<sysN>/{output.dcd,gamd.log,production-start-step.txt,<prmtop>}
# Also keeps config.xml/input.xml/gamd-restart.dat/temperature.dat/checkpoint/runner.log for reproducibility and reweighting.
# Safety: byte-compares every file local vs remote; a mismatch FAILs loudly (no silent pass).
set -uo pipefail
HOST=autodl-nmb1
REMOTE=/root/autodl-tmp/PARPi_design/runs_ruca_prot
DEST="${DATA_ROOT:-/Volumes/tjogzt4T/PARPi_data}/rerun_202609/runs_ruca_prot"
WHICH=${1:?usage: pull_ruca_prot.sh <sys1|sys2|both>}

VOL="${DATA_ROOT:-/Volumes/tjogzt4T}"; [ -d "$VOL" ] || { echo "ABORT: data volume $VOL not mounted"; exit 1; }
mkdir -p "$DEST"

declare -a SRCS
if [ "$WHICH" = "sys1" ] || [ "$WHICH" = "both" ]; then
  S=sys1_rucaparib_prot; SYS=sys1   # archive dir name mirrors runs_nira_prot/sys1
  SRCS+=("$S/gamd_out/output.dcd:$SYS/output.dcd"
         "$S/gamd_out/gamd.log:$SYS/gamd.log"
         "$S/gamd_out/production-start-step.txt:$SYS/production-start-step.txt"
         "$S/gamd_out/gamd-restart.dat:$SYS/gamd-restart.dat"
         "$S/gamd_out/temperature.dat:$SYS/temperature.dat"
         "$S/gamd_out/input.xml:$SYS/input.xml"
         "$S/gamd_out/gamd_restart.checkpoint:$SYS/gamd_restart.checkpoint"
         "$S/config.xml:$SYS/config.xml"
         "$S/runner.log:$SYS/runner.log"
         "$S.prmtop:$SYS/sys1_rucaparib_prot.prmtop")
fi
if [ "$WHICH" = "sys2" ] || [ "$WHICH" = "both" ]; then
  S=sys2_rucaparib_prot; SYS=sys2   # archive dir name mirrors runs_nira_prot/sys2
  SRCS+=("$S/gamd_out/output.dcd:$SYS/output.dcd"
         "$S/gamd_out/gamd.log:$SYS/gamd.log"
         "$S/gamd_out/production-start-step.txt:$SYS/production-start-step.txt"
         "$S/gamd_out/gamd-restart.dat:$SYS/gamd-restart.dat"
         "$S/gamd_out/temperature.dat:$SYS/temperature.dat"
         "$S/gamd_out/input.xml:$SYS/input.xml"
         "$S/gamd_out/gamd_restart.checkpoint:$SYS/gamd_restart.checkpoint"
         "$S/config.xml:$SYS/config.xml"
         "$S/runner.log:$SYS/runner.log"
         "$S/build/rucaparib_best_prot.prmtop:$SYS/rucaparib_best_prot.prmtop")
fi

FAIL=0
for pair in "${SRCS[@]}"; do
  R=${pair%%:*}; L=${pair#*:}
  mkdir -p "$DEST/$(dirname "$L")"
  RS=$(ssh -o ConnectTimeout=15 "$HOST" "stat -c %s $REMOTE/$R" 2>/dev/null)
  [ -z "$RS" ] && { echo "MISS  $R (absent on remote)"; FAIL=1; continue; }
  for try in 1 2 3; do
    rsync -a --size-only --partial -e "ssh -o ConnectTimeout=15" "$HOST:$REMOTE/$R" "$DEST/$L" && break
    echo "  retry $try  $R"
  done
  LS=$(stat -f %z "$DEST/$L" 2>/dev/null || echo 0)
  if [ "$RS" = "$LS" ]; then
    echo "OK    $L  $LS bytes"
  else
    echo "FAIL  $L  local=$LS remote=$RS"; FAIL=1
  fi
done

echo "---- local download directory ----"
du -sh "$DEST"/* 2>/dev/null
echo "---- volume free space ----"; df -h "$VOL" | tail -1
[ "$FAIL" = 0 ] && echo "RESULT: ALL OK" || echo "RESULT: HAS FAILURES"
exit $FAIL
