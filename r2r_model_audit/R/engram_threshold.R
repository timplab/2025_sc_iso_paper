#!/usr/bin/env Rscript
# Reproduce the Fig 4A engram threshold exactly as coded in Figure_4_plotting.Rmd (lines 20-102),
# then compare with a true z > 0.5 threshold (mean + 0.5 SD). Read-only on Sheridan's objects.
# Usage: Rscript engram_threshold.R <out_root>
args <- commandArgs(trailingOnly = TRUE)
out_root <- args[1]
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(DESeq2) })
options(Seurat.object.assay.version = "v3")
DATA_DIR <- "/atium/Data/projects/sheridansinglecell/R2C2/Data"

counts_tx <- readRDS(file.path(DATA_DIR, "CA1/genes_scaled_cond_nobatch.rds"))
combined <- CreateSeuratObject(counts = counts(counts_tx, normalized = FALSE), meta.data = data.frame(colData(counts_tx)))
scalef <- mean(combined@meta.data[["nCount_RNA"]])
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = scalef)
cnt <- GetAssayData(object = combined, slot = "data")
keep <- Matrix::rowSums(cnt) > 0
combined <- CreateSeuratObject(cnt[keep, ], meta.data = combined@meta.data)
fos <- combined@assays[["RNA"]]@counts["Fos", ]          # log-normalized values stored in the counts slot
cells <- colnames(combined)
grp <- ifelse(grepl("UC", cells), "U+", ifelse(grepl("FC", cells), "C+", ifelse(grepl("NC", cells), "N", NA)))
avg <- mean(fos); s <- sd(fos)
thr_coded <- (avg + s) / 2                                 # zpoint5 <- z1/2, z1 <- avg + sd
thr_z05 <- avg + 0.5 * s
eng_coded <- cells[fos > thr_coded]
eng_z05 <- cells[fos > thr_z05]
saved <- readRDS(file.path(DATA_DIR, "CA1/engram_list.rds"))

# Check the raw-count recipe: log1p(Fos / colSum * scalef)
raw <- counts(counts_tx, normalized = FALSE)
fos_manual <- log1p(raw["Fos", ] / colSums(raw) * scalef)

per_group <- function(eng) {
  do.call(rbind, lapply(c("C+", "U+", "N"), function(g) {
    n <- sum(grp == g); k <- sum(cells[grp == g] %in% eng)
    data.frame(group = g, n_cells = n, n_engram = k, pct_engram = 100 * k / n)
  }))
}
res <- list(seurat = as.character(packageVersion("Seurat")), seuratobject = as.character(packageVersion("SeuratObject")),
  deseq2 = as.character(packageVersion("DESeq2")),
  n_cells = length(cells), scalef = scalef, fos_mean = avg, fos_sd = s,
  thr_coded = thr_coded, thr_z05 = thr_z05, z_of_coded_threshold = (thr_coded - avg) / s,
  n_fos_zero = sum(fos == 0), fos_manual_maxdiff = max(abs(fos_manual[cells] - fos)),
  identical_to_saved = identical(sort(eng_coded), sort(saved)),
  n_saved = length(saved), n_coded = length(eng_coded), n_z05 = length(eng_z05),
  n_both = length(intersect(eng_coded, eng_z05)), n_coded_only = length(setdiff(eng_coded, eng_z05)),
  n_z05_only = length(setdiff(eng_z05, eng_coded)),
  by_group_coded = per_group(eng_coded), by_group_z05 = per_group(eng_z05),
  per_cell = data.frame(cell = cells, group = grp, fos_lognorm = as.numeric(fos),
                        engram_coded = cells %in% eng_coded, engram_z05 = cells %in% eng_z05,
                        in_saved_list = cells %in% saved),
  sessionInfo = paste(capture.output(sessionInfo()), collapse = "\n"))
od <- file.path(out_root, "cache", paste0("seurat_", res$seurat)); dir.create(od, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(od, "engram_threshold.rds"))
cat("scalef", scalef, "mean", avg, "sd", s, "thr_coded", thr_coded, "thr_z05", thr_z05,
    "identical_to_saved", res$identical_to_saved, "n_coded", length(eng_coded), "n_z05", length(eng_z05), "\n")
print(res$by_group_coded); print(res$by_group_z05)
