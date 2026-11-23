# Reviewer 3 doublet check (GENOME/2025/281698 revision)

`r2r_doublet_check.Rmd` tests whether the Plp1+ non-engram CA1 cells (Fig 3D, 4C, 4E) are neuron-oligodendrocyte doublets.
It reads Sheridan Cavalier's saved objects read-only and writes only to
`/atium/Data/projects/sheridansinglecell/r2r_doublet_check/`.

Inputs (read-only):
- `/atium/Data/projects/sheridansinglecell/R2C2/Data/Seurat_objs/R2C2_combined_final.rds` (combined ONT object; NC1 is dropped)
- `/atium/Data/projects/sheridansinglecell/R2C2/Data/CA1/genes_scaled_cond_nobatch.rds` (CA1 DESeqDataSet behind Fig 3D)
- `/atium/Data/projects/sheridansinglecell/R2C2/Data/CA1/engram_list.rds` (Fig 4A engram cells)
- Sheridan's doublet lists, for concordance only: `files_for_resub/doublet_removed_data/doublet_list.rds`, plus
  `doublet_list_8perc.rds` and `doublet_list_newpKs.rds` copied unchanged from winry:/home/sherbear/ into `<out>/inputs_copied/`.

Environment (zoro):
```
mamba create -n r2r_doublets -c conda-forge -c bioconda r-base=4.4 r-seurat r-seuratobject bioconductor-deseq2 \
  r-remotes r-rmarkdown r-knitr r-ggplot2 r-patchwork r-cowplot r-ggrepel r-dplyr r-tidyr r-readr r-matrix r-fields r-kernsmooth r-rocr pandoc
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", ref="1b244d8f0d54b4b1cb4365639931bbb16f01e1cd")'
```
Versions used: R 4.4.3, Seurat 5.5.1, SeuratObject 5.4.0, DoubletFinder 2.0.6 (commit 1b244d8), DESeq2 1.46.0 (full sessionInfo in the output dir).

Run (on a host that mounts /atium):
```
export PATH=~/miniforge3/envs/r2r_doublets/bin:$PATH
Rscript -e 'rmarkdown::render("r2r_doublet_check.Rmd", output_dir="/atium/Data/projects/sheridansinglecell/r2r_doublet_check", params=list(ncores=32))'
```
DoubletFinder calls and DESeq2 refits are cached in `<out>/cache/`; pass `params=list(force=TRUE)` to recompute.
First run about 30 min on zoro; cached reruns about 3 min.

Outputs: `FigS_doublet_check.pdf/.png`, `source_data_figS/` (one CSV per panel), `tables/`, `r2r_doublet_check.html`, `sessionInfo.txt`.
(`source_data/` in the output dir is from a superseded draft with different panel letters; use `source_data_figS/`.)
