# R2R model audit (GENOME/2025/281698)

Audit of the DESeq2 models behind the experience-variable gene and isoform lists (Supp Data 6, Figs 3 and 6), the engram
models (Fig 4C/4D/4E, Fig S21B/C), the cell-type model (Fig 5, Supp Data 4) and the time models (Fig S9, Reviewer 2 #3),
plus clean refits. Everything reads Sheridan Cavalier's saved objects **read-only** from
`/atium/Data/projects/sheridansinglecell/R2C2/Data/` and writes only to
`/atium/Data/projects/sheridansinglecell/r2r_model_audit/`.

## Layout

- `R/audit_functions.R` shared helpers: clean refit, Sheridan's wrap-and-refit call, results extraction, experience labels,
  size-factor provenance, filter statistics.
- `R/run_job.R` one job = one model on one saved object; writes `cache/<DESeq2 version>/<analysis>__<object>.rds`.
  `job_table()` lists every job; `objects` maps object keys to files.
- `R/engram_threshold.R` reproduces the Fig 4A engram threshold exactly as coded (Figure_4_plotting.Rmd) and the z > 0.5 alternative.
- `run_fits.sh` runs a job list on one host with `xargs -P`. Finished jobs are skipped (cached).
- `jobs/` job lists used for the run (`w_*` winry/1.38.3, `z_*` zoro/1.46.0, each split into small and big objects).
- `r2r_model_audit.Rmd` assembles the cached fits into the tables (`tables/table1` to `table7`), source data and sessionInfo.
  It fits nothing itself.
- `render.sh` renders the Rmd on zoro.
- `launch_zoro.sh`, `rerun_oom.sh` the exact launch and out-of-memory rerun wrappers used for this run.
- `deliver/` the one-page summary (Case A/B verdict) and draft Methods, legend, Results and reply text.

## Software

| role | host | R | DESeq2 | Seurat |
|---|---|---|---|---|
| reproduce Sheridan (her library, read-only via `R_LIBS_USER`) | winry | 4.2.2 | 1.38.3 | 5.3.0 |
| current release | zoro, conda env `r2r_doublets` | 4.4.3 | 1.46.0 | 5.5.1 |

/atium is mounted on zoro and winry only. Sheridan's library lives on winry
(`/home/sherbear/R/x86_64-pc-linux-gnu-library/4.2`), so all 1.38.3 jobs run there.

## How to run

```bash
OUT=/atium/Data/projects/sheridansinglecell/r2r_model_audit
mkdir -p $OUT && rsync -a r2r_model_audit/ $OUT/code/
# 1. on winry (DESeq2 1.38.3, Sheridan's library)
cd $OUT/code
./run_fits.sh 1.38.3 jobs/w_1.38.3_big.txt 2 10 &
./run_fits.sh 1.38.3 jobs/w_1.38.3_small.txt 6 8 &
R_LIBS_USER=/home/sherbear/R/x86_64-pc-linux-gnu-library/4.2 /usr/bin/Rscript R/engram_threshold.R $OUT
wait
# 2. on zoro (DESeq2 1.46.0)
cd $OUT/code
./run_fits.sh 1.46.0 jobs/z_1.46.0_big.txt 3 8 &
./run_fits.sh 1.46.0 jobs/z_1.46.0_small.txt 6 8 &
wait
# 3. on zoro: assemble tables and HTML
./render.sh
```

Arguments to `run_fits.sh`: DESeq2 version, job list, number of jobs at once, cores per job (DESeq `parallel=TRUE`,
`MulticoreParam`). The first full run took about 2.5 h on winry and 1.3 h on zoro (winry ran out of memory once with five large
isoform jobs at a time; keep the big list at 2 jobs x 10 cores there). Re-running skips cached jobs; delete a
cache file to refit that one model.

## Models fitted (all LRT, `fitType="local"`, `minmu=1e-6`, `minReplicatesForReplace=Inf`, `useT=TRUE` as in Sheridan's calls)

| analysis | full | reduced | notes |
|---|---|---|---|
| `stored` | (none) | | `results()` on the saved fit, plus the inventory |
| `exp_sheridan` | `DESeqDataSet(saved, ~ sex + batch + condition)` then `DESeq()` | `~ sex + batch` | Sheridan's call as written |
| `exp_scaled` | `~ scaled + condition` | `~ scaled` | clean refit; reproduces the saved fit |
| `exp_sbc` | `~ sex + batch + condition` | `~ sex + batch` | clean refit of the model the Methods describe |
| `exp_ssbc` | `~ scaled + sex + batch + condition` | `~ scaled + sex + batch` | added sensitivity model (keeps the library-size covariate) |
| `exp_*_N` | as above | | Naive as baseline instead of Conditioned (CA1 only) |
| `engram_sheridan` | wrap of saved CA1 fit, `~ sex + batch + engram` | `~ sex + batch` | Figure_4 / Figure_7 as written |
| `engram_sbe` | `~ sex + batch + engram` | `~ sex + batch` | clean |
| `engram_sbce` | `~ sex + batch + condition + engram` | `~ sex + batch + condition` | clean, adjusts for experience |
| `timeC_*`, `timeU_*` | `~ sex + time` within C+ or U+ cells | `~ sex` | plotting_for_revision.Rmd; `_sheridan` as written, `_clean` refit |
| `condtime_*` | `~ sex + batch + condition_time` (Con0, Con1, Ucon0, Ucon1, Naive) | `~ sex + batch` | condition_time_for_supp.Rmd (winry:/home/sherbear) |
| `celltype_bt` | `~ batch + type` on naive cells | `~ batch` | clean refit of the Fig 5 cell-type model; `~ sex + batch + type` is not full rank |

"Clean" means a new `DESeqDataSetFromMatrix()` built from the saved counts and colData, with the saved size factors.
Saved size factors are 1 for every isoform object; for gene objects they are DESeq2 median-of-ratios (`type="ratio"`)
computed on the 1 to 5 genes detected in every cell of that cell type.
