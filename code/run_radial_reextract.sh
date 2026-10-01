#!/bin/bash
# Re-extract the radial representation: benchmark first, then the catalogue.
#
# Separate from run_radial.sh, which drives the original stage-2 job over its
# own worklist. This one re-extracts everything through the post-2026-10-01
# fetch path.
#
# The order is the point. The SIA service stopped serving 0.05 deg boxes, and
# re-extracting objects already in radial_profiles2.npz reproduced them only to
# corr 0.963 with 18.5% median relative difference. So the model trained on the
# August profiles cannot be trusted against anything extracted now. Doing the
# 24,049 benchmark objects FIRST means the representation can be retrained and
# checked after about two days; if the AUC does not survive the new extraction,
# that is known before another thirteen days go into the catalogue.
#
# Both stages are shard-resumable. Kill this freely and restart; finished
# shards are skipped.
#
#   setsid nohup ./sweeps/run_radial_reextract.sh \
#       > /home/ben/legacy/radial_reextract.log 2>&1 &
#
# Progress:
#   .venv/bin/python sweeps/radial_catalogue.py --benchmark --status
#   .venv/bin/python sweeps/radial_catalogue.py --status
set -u

# Same lock discipline as the other wrappers, with its own lock file so this
# cannot collide with run_radial.sh.
exec 9>/home/ben/legacy/.radial_reextract.lock
flock -n 9 || { echo "already running"; exit 0; }
echo $$ > /home/ben/legacy/radial_reextract.pid

cd /home/ben/star
PY=/home/ben/star/.venv/bin/python

echo "=== $(date -Is)  stage 1 of 2: benchmark, 24,049 objects"
"$PY" -u sweeps/radial_catalogue.py --benchmark
echo "=== $(date -Is)  benchmark stage finished"
"$PY" -u sweeps/radial_catalogue.py --benchmark --status

# Stop unless the benchmark stage actually produced something. Running the
# catalogue on an extraction path that yielded nothing would waste a fortnight.
BENCH=$(ls /home/ben/legacy/radial_bench/shard_*.npz 2>/dev/null | wc -l)
if [ "$BENCH" -lt 10 ]; then
  echo "=== only ${BENCH}/13 benchmark shards written; NOT starting the catalogue."
  echo "=== investigate, then rerun this script to resume."
  exit 1
fi

echo "=== $(date -Is)  stage 2 of 2: catalogue, 161,721 objects"
"$PY" -u sweeps/radial_catalogue.py
echo "=== $(date -Is)  catalogue stage finished"
"$PY" -u sweeps/radial_catalogue.py --status
