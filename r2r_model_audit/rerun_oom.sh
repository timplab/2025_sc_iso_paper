#!/usr/bin/env bash
cd /atium/Data/projects/sheridansinglecell/r2r_model_audit/code
while pgrep -f '^/usr/lib/R/bin/exec/R.*run_job.R' > /dev/null || pgrep -f 'run_job.R exp_ssbc|run_job.R timeC' > /dev/null; do sleep 30; done
printf "exp_ssbc__txgc_Dentate_1\ncelltype_bt__allnaive_type\n" > jobs/w_rerun_oom.txt
./run_fits.sh 1.38.3 jobs/w_rerun_oom.txt 1 10

