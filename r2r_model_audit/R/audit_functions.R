# Shared functions for the R2R DE model audit (GENOME/2025/281698).
# Works under DESeq2 1.38.3 (Sheridan's library, R 4.2.2, winry) and 1.46.0 (r2r_doublets env, R 4.4.3, zoro).
# Reads Sheridan's saved objects read-only. Never writes to her directories.

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
})

DATA_DIR <- "/atium/Data/projects/sheridansinglecell/R2C2/Data"

# DESeq() settings used in every one of Sheridan's calls. useT only affects Wald tests (DESeq() does not
# pass it to nbinomLRT), so it is inert for these LRTs; we keep it so the calls match hers.
deseq_lrt <- function(dds, reduced, ncores = 1) {
  if (ncores > 1) {
    DESeq(dds, fitType = "local", useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf,
          test = "LRT", reduced = reduced, parallel = TRUE, BPPARAM = MulticoreParam(ncores))
  } else {
    DESeq(dds, fitType = "local", useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf,
          test = "LRT", reduced = reduced)
  }
}

# Clean object: counts + colData (minus sizeFactor) of the saved object, fresh rowData, stored size factors.
clean_from <- function(d, design, relevel_naive = FALSE) {
  cd <- colData(d)
  cd$sizeFactor <- NULL
  if (relevel_naive) cd$condition <- relevel(factor(cd$condition), "Naive")
  dd <- suppressWarnings(suppressMessages(
    DESeqDataSetFromMatrix(countData = counts(d), colData = cd, design = design)))
  sizeFactors(dd) <- sizeFactors(d)
  dd
}

# Sheridan's pattern, as written: wrap the (already fitted) object and call DESeq() again.
sheridan_wrap <- function(d, design) {
  suppressWarnings(suppressMessages(DESeqDataSet(d, design = design)))
}

# One named mcols column; when the name is duplicated return the first and last copies.
mcol_first_last <- function(dds, name) {
  mc <- mcols(dds)
  idx <- which(colnames(mc) == name)
  if (length(idx) == 0) return(list(first = rep(NA_real_, nrow(dds)), last = rep(NA_real_, nrow(dds)), n = 0L))
  list(first = mc[[idx[1]]], last = mc[[idx[length(idx)]]], n = length(idx))
}

# Cook's cutoff that results() uses by default: qf(0.99, p, m - p).
cooks_cutoff_default <- function(dds) {
  p <- ncol(attr(dds, "modelMatrix"))
  m <- ncol(dds)
  qf(0.99, p, m - p)
}

# Per-feature table. contrasts: named list of c(factor, numerator, denominator).
extract_results <- function(dds, contrasts = list()) {
  r <- results(dds)
  az <- mcols(dds)$allZero
  if (is.null(az)) az <- rowSums(counts(dds)) == 0
  mc <- mcol_first_last(dds, "maxCooks")
  pv <- mcol_first_last(dds, "LRTPvalue")
  out <- data.frame(feature = rownames(dds),
                    baseMean = r$baseMean,
                    stat = r$stat,
                    pvalue = r$pvalue,
                    padj = r$padj,
                    allZero = az,
                    maxCooks_first = mc$first,
                    maxCooks_last = mc$last,
                    LRTPvalue_first = pv$first,
                    LRTPvalue_last = pv$last,
                    stringsAsFactors = FALSE)
  for (nm in names(contrasts)) {
    rc <- results(dds, contrast = contrasts[[nm]])
    out[[paste0("lfc_", nm)]] <- rc$log2FoldChange
    attr(out, paste0("padj_equal_", nm)) <- isTRUE(all.equal(rc$padj, r$padj))
  }
  out
}

fit_meta <- function(dds, label) {
  mc <- mcols(dds)
  mm <- mcols(mc)
  dups <- unique(colnames(mc)[duplicated(colnames(mc))])
  lrt_desc <- mm$description[colnames(mc) == "LRTPvalue"]
  list(label = label,
       deseq2_version = as.character(packageVersion("DESeq2")),
       design = paste(deparse(design(dds)), collapse = ""),
       lrt_description = paste(lrt_desc, collapse = " || "),
       n_features = nrow(dds), n_cells = ncol(dds),
       dup_mcols = paste(dups, collapse = ";"),
       n_dup_mcols = length(dups),
       cooks_cutoff = tryCatch(cooks_cutoff_default(dds), error = function(e) NA_real_),
       model_matrix_cols = paste(colnames(attr(dds, "modelMatrix")), collapse = ";"),
       reduced_matrix_cols = paste(colnames(attr(dds, "reducedModelMatrix")), collapse = ";"),
       n_fullBetaConv_false = sum(!mcols(dds)$fullBetaConv, na.rm = TRUE),
       fitType = attr(dispersionFunction(dds), "fitType"))
}

# Sheridan's experience assignment (DESeq2_analysis.Rmd), vectorised. Needs padj and the three
# pairwise LFCs (C vs U, C vs N, U vs N); the other three directions are their negatives.
assign_experience <- function(padj, lfc_CvU, lfc_CvN, lfc_UvN, lfc_thr = 0.25, q_thr = 0.05) {
  gt <- function(x) !is.na(x) & x > lfc_thr
  sig <- !is.na(padj) & padj < q_thr
  cond_up   <- sig & (gt(lfc_CvU)  | gt(lfc_CvN))
  uncond_up <- sig & (gt(-lfc_CvU) | gt(lfc_UvN))
  naive_up  <- sig & (gt(-lfc_CvN) | gt(-lfc_UvN))
  lab <- rep(NA_character_, length(padj))
  lab[sig] <- "Log2 < .25"
  lab[cond_up & !uncond_up & !naive_up] <- "C+"
  lab[uncond_up & !cond_up & !naive_up] <- "U+"
  lab[naive_up & !cond_up & !uncond_up] <- "N"
  lab[cond_up & uncond_up] <- "C+U+"
  lab[cond_up & naive_up] <- "C+N"
  lab[uncond_up & naive_up] <- "U+N"
  lab
}

# Size-factor provenance: which simple rule reproduces the stored size factors (up to a constant)?
sf_provenance <- function(d) {
  sf <- sizeFactors(d)
  if (is.null(sf)) return(list(sf_stored = FALSE))
  gm <- function(x) exp(mean(log(x)))
  cmp <- function(a) { if (any(!is.finite(a)) || any(a <= 0)) return(NA_real_); max(abs(log(sf / gm(sf)) - log(a / gm(a)))) }
  cnt <- counts(d)
  cs <- colSums(cnt)
  n_allpos <- sum(rowSums(cnt > 0) == ncol(cnt))
  ratio <- tryCatch(suppressWarnings(estimateSizeFactorsForMatrix(cnt)), error = function(e) rep(NA_real_, ncol(cnt)))
  pos <- tryCatch(estimateSizeFactorsForMatrix(cnt, type = "poscounts"), error = function(e) rep(NA_real_, ncol(cnt)))
  ncr <- if (!is.null(d$nCount_RNA)) d$nCount_RNA else rep(NA_real_, ncol(cnt))
  list(sf_stored = TRUE, sf_all_one = all(sf == 1),
       sf_min = min(sf), sf_median = median(sf), sf_max = max(sf),
       n_features_positive_in_all_cells = n_allpos,
       sf_maxlogdiff_ratio = cmp(ratio), sf_maxlogdiff_poscounts = cmp(pos),
       sf_maxlogdiff_colsums = cmp(cs), sf_maxlogdiff_nCount_RNA = cmp(ncr))
}

# Detection-based filter stats for one count matrix.
filter_stats <- function(cnt) {
  ncell <- rowSums(cnt > 0)
  tot <- rowSums(cnt)
  nz <- ncell > 0
  list(n_features_matrix = nrow(cnt), n_allzero = sum(!nz), n_detected_any = sum(nz),
       min_cells_detected_nonzero = if (any(nz)) min(ncell[nz]) else NA,
       min_total_count_nonzero = if (any(nz)) min(tot[nz]) else NA,
       n_detected_ge3 = sum(ncell >= 3), n_detected_ge10 = sum(ncell >= 10),
       q10_cells_detected_nonzero = if (any(nz)) unname(quantile(ncell[nz], 0.1)) else NA,
       median_cells_detected_nonzero = if (any(nz)) median(ncell[nz]) else NA)
}
