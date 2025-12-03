#!/usr/bin/env bash

# Number of parallel jobs
NJOBS=30

# Name of the ROOT macro file (must contain int SigmaProton())
MACRO=SigmaProton.C

# Base output directory
OUTBASE="SigmaProton_runs"

mkdir -p "${OUTBASE}"

for i in $(seq 0 $((NJOBS - 1))); do
  outdir="${OUTBASE}/run_${i}"
  mkdir -p "${outdir}"

  # Copy the macro into the job directory (or you can just reference it by path)
  cp "${MACRO}" "${outdir}/"

  (
    cd "${outdir}"

    # Optional: set a different seed per job via env var (if you adapt the macro to use it)
    export JOBID="${i}"

    # Run the macro in batch mode; ACLiC will compile it the first time
    # If you want to force recompilation, use SigmaProton.C++
    root -l -b -q 'SigmaProton.C+' &> job.log
  ) &
done

# Wait for all background jobs to finish
wait

echo "All ${NJOBS} jobs finished. Outputs are in ${OUTBASE}/run_*/sigmaprotontree_test.root"
