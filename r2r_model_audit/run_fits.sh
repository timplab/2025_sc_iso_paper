#!/usr/bin/env bash
# Run a list of audit jobs on one host.
# usage: run_fits.sh <1.38.3|1.46.0> <job_list.txt> <n_parallel_jobs> <cores_per_job>
#   1.38.3: run on winry (Sheridan's R 4.2.2 library, read-only via R_LIBS_USER)
#   1.46.0: run on zoro  (conda env r2r_doublets, R 4.4.3)
set -u
VER=$1; JOBS=$2; NP=$3; NC=$4
OUT=${AUDIT_OUT:-/atium/Data/projects/sheridansinglecell/r2r_model_audit}
HERE=$(cd "$(dirname "$0")" && pwd)
if [ "$VER" = "1.38.3" ]; then
  export R_LIBS_USER=/home/sherbear/R/x86_64-pc-linux-gnu-library/4.2; RS=/usr/bin/Rscript
else
  export PATH=$HOME/miniforge3/envs/r2r_doublets/bin:$PATH; RS=Rscript
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$OUT/logs/$VER"
grep -v '^#' "$JOBS" | grep -v '^$' | xargs -P "$NP" -I{} sh -c "$RS $HERE/R/run_job.R {} $OUT $NC > $OUT/logs/$VER/{}.log 2>&1; echo \"\$(date +%F_%T) {} exit \$?\" >> $OUT/logs/$VER/_progress.txt"
