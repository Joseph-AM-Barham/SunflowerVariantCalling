#!/usr/bin/env bash
# SLURM job status checker for Snakemake --cluster-status.
# Prints "running", "success", or "failed" for the given job ID.
# Called by Snakemake after each cluster submission to track job state.

jobid="$1"

# sacct may not be populated immediately after submission; fall back to squeue.
state=$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null \
        | head -1 | tr -d ' ')

if [ -z "$state" ]; then
    # sacct not yet populated -- check squeue instead
    state=$(squeue -j "$jobid" -h -o "%t" 2>/dev/null)
    case "$state" in
        R|PD|CG) echo "running" ;;
        "")      echo "success" ;;  # not in queue and not in sacct -- assume done
        *)       echo "running" ;;
    esac
    exit 0
fi

case "$state" in
    COMPLETED)                              echo "success" ;;
    FAILED|CANCELLED*|TIMEOUT|NODE_FAIL)    echo "failed"  ;;
    RUNNING|PENDING|REQUEUED|COMPLETING)    echo "running" ;;
    "")                                     echo "running" ;;
    *)                                      echo "running" ;;
esac
