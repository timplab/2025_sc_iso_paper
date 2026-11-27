#!/usr/bin/env bash
cd /atium/Data/projects/sheridansinglecell/r2r_model_audit/code
./run_fits.sh 1.46.0 jobs/z_1.46.0_big.txt 3 8 > /atium/Data/projects/sheridansinglecell/r2r_model_audit/logs/z_big.out 2>&1 &
./run_fits.sh 1.46.0 jobs/z_1.46.0_small.txt 6 8 > /atium/Data/projects/sheridansinglecell/r2r_model_audit/logs/z_small.out 2>&1 &
(export PATH=$HOME/miniforge3/envs/r2r_doublets/bin:$PATH; Rscript R/engram_threshold.R /atium/Data/projects/sheridansinglecell/r2r_model_audit > /atium/Data/projects/sheridansinglecell/r2r_model_audit/logs/engram_threshold_zoro.out 2>&1) &
wait

