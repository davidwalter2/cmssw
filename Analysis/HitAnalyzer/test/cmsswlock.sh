#!/bin/bash
# Guard against `scram b` replacing a plugin under a running cmsRun.
#
# WHY (2026-08-15). Two cmsRun jobs died during the profiling work with what
# looked like physics crashes -- one SIGSEGV inside
# RunManagerMTWorker::beginRun, one SIGBUS mid-event-loop. Neither was:
# biglib/el9_amd64_gcc12/pluginSimulation.so was rewritten at 10:29:22 and
# again at 10:33:03 by a concurrent `scram b` in the same shared area, and a
# process whose mmap'd shared object is replaced underneath it takes SIGBUS on
# the next page fault into it. That presents as a mysterious segfault with a
# corrupted stack (`#4 0x00000000000f1070 in ?? ()`), which is exactly how much
# time it costs to diagnose.
#
# This is an advisory lock, so it only works if EVERYONE uses it -- hence the
# two modes:
#
#   cmsswlock.sh run   <cmd...>    SHARED    -- any number of concurrent jobs
#   cmsswlock.sh build <cmd...>    EXCLUSIVE -- waits for running jobs to finish
#
# A build therefore waits for the current jobs, and a job started once the build
# HOLDS the lock waits for the build. Neither ever sees a half-written .so.
#
# Measured behaviour (2026-08-15): two `run` holders overlapped; a `build`
# requested while they held it started 4 ms after the last one finished. The one
# gap is flock's: a `run` requested while a `build` is merely QUEUED is still
# granted immediately, so a continuous stream of jobs can starve a build. That
# delays the build, it does not corrupt anything -- if you are being starved,
# stop submitting or use a private area.
#
# The alternative, and the better answer for long campaigns, is a PRIVATE build
# area per worker -- the convention the pT-scan work used with `toyscan_pt/`.
# Use that when a campaign will run for hours; use this lock for the ordinary
# case of "someone else is also working in this area right now".
#
# usage:
#   ./cmsswlock.sh run   cmsRun runToyGeomCheck.py events=2000
#   ./cmsswlock.sh build scram b -j 32
#   LOCK_TIMEOUT=1800 ./cmsswlock.sh build scram b        # default 3600 s
set -uo pipefail

MODE=${1:-}
shift || true
if [[ -z "$MODE" || $# -eq 0 ]]; then
    echo "usage: $0 {run|build} <command...>" >&2
    exit 2
fi

AREA=${CMSSW_BASE:-}
if [[ -z "$AREA" ]]; then
    echo "$0: CMSSW_BASE is not set -- run cmsenv first" >&2
    exit 2
fi
LOCK="$AREA/.plugin_build.lock"
: > "$LOCK" 2>/dev/null || true
TIMEOUT=${LOCK_TIMEOUT:-3600}

case "$MODE" in
    run)   FLAG=-s; WHAT="shared (runner)" ;;
    build) FLAG=-x; WHAT="EXCLUSIVE (build)" ;;
    *)     echo "$0: mode must be 'run' or 'build', got '$MODE'" >&2; exit 2 ;;
esac

exec 9<>"$LOCK"
if ! flock "$FLAG" -w "$TIMEOUT" 9; then
    echo "$0: timed out after ${TIMEOUT}s waiting for the $WHAT lock on $LOCK" >&2
    echo "     (something is holding it: check with 'fuser -v $LOCK')" >&2
    exit 75
fi
"$@"
rc=$?
exec 9>&-
exit $rc
