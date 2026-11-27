#!/usr/bin/env bash
export PATH=$HOME/miniforge3/envs/r2r_doublets/bin:$PATH
OUT=/atium/Data/projects/sheridansinglecell/r2r_model_audit
mkdir -p $OUT/render_tmp
Rscript -e "rmarkdown::render('$OUT/code/r2r_model_audit.Rmd', output_dir='$OUT', intermediates_dir='$OUT/render_tmp', knit_root_dir='$OUT/render_tmp')" > $OUT/logs/render.log 2>&1
echo "render exit $?" >> $OUT/logs/render.log
