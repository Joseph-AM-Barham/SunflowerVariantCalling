#!/usr/bin/env bash
# PBS job status checker for Snakemake --cluster-status.
# Prints "running", "success", or "failed" for the given job ID.
# Called by Snakemake after each cluster submission to track job state.

jobid="$1"

state=$(qstat -f "$jobid" 2>/dev/null | awk '/job_state[[:space:]]*=/ {print $3}')

case "$state" in
    R|Q|H|T|W|S) echo "running"  ;;   # Running, Queued, Held, Moving, Waiting, Suspended
    C)            echo "success"  ;;   # Completed
    E|F)          echo "failed"   ;;   # Exiting (error), Failed
    "")           echo "success"  ;;   # No longer in queue -- assume completed successfully
    *)            echo "running"  ;;   # Unknown state -- assume still running
esac
