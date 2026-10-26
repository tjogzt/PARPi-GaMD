#!/bin/bash
# watch_tala_26ns.sh — local watchdog for the tala S2 26-ns run.
# Polls every 15 min; exits 0 when the DCD reaches 760 frames (38 ns total,
# 26 ns production), exits 1 if the runner disappears before completion.
set -u
PWFILE=/Users/taozhu/.hermes/cache/scratch/.autodl_pw
HOST=connect.weste.seetacloud.com
PORT=47874
RDIR=/root/autodl-tmp/PARPi_design/runs_extension/tala_s2_26ns
LOG=/Users/taozhu/clacky_workspace/PARPi_design/review-stage/tala_26ns/watch.log
mkdir -p "$(dirname "$LOG")"
echo "[$(date '+%F %T')] watchdog start" >> "$LOG"
for i in $(seq 1 96); do
  OUT=$(sshpass -f "$PWFILE" ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    -o ConnectTimeout=20 -p "$PORT" root@"$HOST" \
    "cd $RDIR || exit 9
     SIZE=\$(stat -c%s gamd_out/output.dcd 2>/dev/null || echo 0)
     if [ \"\$SIZE\" -gt 92 ]; then FR=\$(( (SIZE-92)/(12*285796+8) )); else FR=0; fi
     STEP=\$(awk '{for(i=1;i<=NF;i++) if(\$i ~ /^[0-9]+\$/ && length(\$i)>=8){print \$i; exit}}' gamd_out/gamd.log 2>/dev/null)
     PROC=\$(pgrep -f '[c]onfig_tala_s2_26ns_seed2026.xml' | head -1)
     echo \"\$FR \${STEP:-0} \${PROC:-none}\"" 2>/dev/null | tail -1)
  RC=$?
  if [ $RC -ne 0 ] || [ -z "${OUT:-}" ]; then
    echo "[$(date '+%F %T')] poll $i: SSH/parse failed" >> "$LOG"
    sleep 900; continue
  fi
  FR=$(echo "$OUT" | awk '{print $1}')
  STEP=$(echo "$OUT" | awk '{print $2}')
  PROC=$(echo "$OUT" | awk '{print $3}')
  echo "[$(date '+%F %T')] poll $i: frames=$FR step=$STEP proc=$PROC" >> "$LOG"
  if [ "${FR:-0}" -ge 760 ]; then
    echo "[$(date '+%F %T')] DONE frames=$FR" >> "$LOG"
    exit 0
  fi
  if [ "$PROC" = "none" ] && [ "${FR:-0}" -lt 760 ]; then
    echo "[$(date '+%F %T')] FAILED: runner gone at frames=$FR" >> "$LOG"
    exit 1
  fi
  sleep 900
done
echo "[$(date '+%F %T')] TIMEOUT after 24h" >> "$LOG"
exit 2
