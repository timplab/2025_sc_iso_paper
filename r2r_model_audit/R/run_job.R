#!/usr/bin/env Rscript
# Run one audit job and cache a compact result.
# Usage: Rscript run_job.R <job_name> <out_root> [ncores]
#   writes <out_root>/cache/<DESeq2 version>/<job_name>.rds
# Job names are <analysis>__<object_key>; see job_table().

args <- commandArgs(trailingOnly = TRUE)
script_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else Sys.getenv("AUDIT_R_DIR", "R")
})
source(file.path(script_dir, "audit_functions.R"))

# ---------------------------------------------------------------- objects
celltypes_gene <- c(CA1 = "CA1/genes_scaled_cond_nobatch.rds", CA2 = "CA2/genes_scaled_cond_nobatch.rds",
  CA3 = "CA3/genes_scaled_cond_nobatch.rds", Astrocyte = "Astrocyte/genes_scaled_cond_nobatch.rds",
  Dentate_1 = "Dentate_1/genes_scaled_cond_nobatch.rds", Dentate_2 = "Dentate_2/genes_scaled_cond_nobatch.rds",
  Dentate_3 = "Dentate_3/genes_scaled_cond_nobatch.rds", Endothelial = "Endothelial/genes_scaled_cond_nobatch.rds",
  Inhibitory = "Inhibitory/genes_scaled_cond_nobatch.rds", Microglia = "Microglia_1/genes_scaled_cond_nobatch.rds",
  Neuroblast = "Neuroblast/genes_scaled_cond_nobatch.rds", Olig = "Olig_scaled_cond_nobatch.rds",
  OPC = "OPC/genes_scaled_cond_nobatch.rds", Pericyte = "Pericyte/genes_scaled_cond_nobatch.rds",
  Progenitor = "Progenitor/genes_scaled_cond_nobatch.rds")
txsc_dirs <- c(CA1 = "CA1", CA2 = "CA2", CA3 = "CA3", Astrocyte = "Astrocyte", Dentate_1 = "Dentate_1",
  Dentate_2 = "Dentate_2", Dentate_3 = "Dentate_3", Endothelial = "Endothelial", Inhibitory = "Inhibitory",
  Microglia = "Microglia_1", Neuroblast = "Neuroblast", Olig_1 = "Olig_1", OPC = "OPC", Pericyte = "Pericyte",
  Progenitor = "Progenitor")
txgc_dirs <- c(txsc_dirs, Olig_2 = "Olig_2")

objects <- rbind(
  data.frame(key = paste0("gene_", names(celltypes_gene)), celltype = names(celltypes_gene), level = "gene",
             file = unname(celltypes_gene)),
  data.frame(key = paste0("txsc_", names(txsc_dirs)), celltype = names(txsc_dirs), level = "isoform",
             file = paste0(txsc_dirs, "/tx_scaled_cond_nobatch.rds")),
  data.frame(key = paste0("txgc_", names(txgc_dirs)), celltype = names(txgc_dirs), level = "isoform",
             file = paste0(txgc_dirs, "/tx_se_gencode_LRT.rds")),
  data.frame(key = "allnaive_type", celltype = "all_naive", level = "isoform",
             file = "allnaivecells_isoform_by_type.rds"))
objects$path <- file.path(DATA_DIR, objects$file)

# isoform objects used for the published experience lists (DESeq2_analysis.Rmd) and the five-level time model
# (condition_time_for_supp.Rmd): tx_scaled_cond_nobatch for CA1 and Inhibitory, tx_se_gencode_LRT otherwise.
exp_iso_keys <- c("txsc_CA1", "txgc_CA2", "txgc_CA3", "txgc_Astrocyte", "txgc_Dentate_1", "txgc_Dentate_2",
  "txgc_Dentate_3", "txgc_Endothelial", "txgc_Microglia", "txgc_Neuroblast", "txgc_Olig_1", "txgc_Olig_2",
  "txsc_Inhibitory", "txgc_OPC", "txgc_Pericyte", "txgc_Progenitor")
exp_gene_keys <- paste0("gene_", names(celltypes_gene))
exp_keys <- c(exp_gene_keys, exp_iso_keys)
# isoform objects used by the within-condition time models (plotting_for_revision.Rmd)
time_iso_keys <- c(paste0("txsc_", names(txsc_dirs)), "txgc_Olig_2")
time_keys <- c(exp_gene_keys, time_iso_keys)

job_table <- function() {
  J <- function(analysis, keys) data.frame(job = paste0(analysis, "__", keys), analysis = analysis, key = keys)
  rbind(
    J("stored", objects$key),
    J("exp_sheridan", exp_keys),
    J("exp_scaled", exp_keys), J("exp_sbc", exp_keys), J("exp_ssbc", exp_keys),
    J("exp_scaled_N", c("gene_CA1", "txsc_CA1")), J("exp_sbc_N", c("gene_CA1", "txsc_CA1")),
    J("engram_sheridan", c("gene_CA1", "txsc_CA1")),
    J("engram_sbe", c("gene_CA1", "txsc_CA1")), J("engram_sbce", c("gene_CA1", "txsc_CA1")),
    J("timeC_sheridan", time_keys), J("timeU_sheridan", time_keys),
    J("timeC_clean", time_keys), J("timeU_clean", time_keys),
    J("condtime_sheridan", exp_keys), J("condtime_clean", exp_keys),
    J("celltype_bt", "allnaive_type"))
}

exp_contrasts <- list(CvU = c("condition", "Conditioned", "Unconditioned"),
                      CvN = c("condition", "Conditioned", "Naive"),
                      UvN = c("condition", "Unconditioned", "Naive"))

engram_cells_for <- function(d) {
  e <- readRDS(file.path(DATA_DIR, "CA1/engram_list.rds"))
  if (!all(e %in% colnames(d))) {
    # isoform objects use FC1_<bc>; gene objects use FC1_FC1_<bc> (Figure_7_plotting.Rmd format())
    e <- vapply(strsplit(e, "_"), function(x) paste0(x[2], "_", x[3]), "")
  }
  stopifnot(all(e %in% colnames(d)))
  e
}

design_rank <- function(cd, f) {
  mm <- tryCatch(stats::model.matrix(f, data = as.data.frame(cd)), error = function(e) e)
  if (inherits(mm, "error")) return(list(ncol = NA, rank = NA, aliased = paste("model.matrix error:", conditionMessage(mm))))
  q <- qr(mm)
  aliased <- if (q$rank < ncol(mm)) colnames(mm)[q$pivot[(q$rank + 1):ncol(mm)]] else character(0)
  list(ncol = ncol(mm), rank = q$rank, aliased = paste(aliased, collapse = ";"))
}

# ---------------------------------------------------------------- run
if (sys.nframe() == 0L && length(args) >= 2) {
  job <- args[1]; out_root <- args[2]; ncores <- if (length(args) >= 3) as.integer(args[3]) else 1L
  jt <- job_table(); stopifnot(job %in% jt$job)
  analysis <- jt$analysis[jt$job == job]; key <- jt$key[jt$job == job]
  ob <- objects[objects$key == key, ]
  ver <- as.character(packageVersion("DESeq2"))
  out_dir <- file.path(out_root, "cache", ver); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(job, ".rds"))
  if (file.exists(out_file)) { message("cached: ", out_file); quit(save = "no") }
  t0 <- Sys.time()
  message(format(t0), " ", job, " DESeq2 ", ver, " ncores ", ncores)
  d <- readRDS(ob$path)
  res <- list(job = job, analysis = analysis, key = key, celltype = ob$celltype, level = ob$level,
              object = ob$path, deseq2 = ver, warnings = character(0))
  wh <- function(expr) withCallingHandlers(expr, warning = function(w) {
    res$warnings <<- c(res$warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  isexp <- "condition" %in% all.vars(design(d))   # allnaive has a constant condition column but a cell-type design

  if (analysis == "stored") {
    mc <- mcols(d); mm <- mcols(mc)
    lrt_desc <- mm$description[colnames(mc) == "LRTPvalue"][1]
    cd <- as.data.frame(colData(d))
    grp <- if (isexp) "condition" else tail(all.vars(design(d)), 1)
    inv <- c(list(object = ob$path, file_bytes = file.size(ob$path), class = class(d)[1],
      n_features = nrow(d), n_cells = ncol(d),
      fitted = all(c("dispersion", "LRTPvalue") %in% colnames(mc)) || "WaldPvalue" %in% colnames(mc),
      design = paste(deparse(design(d)), collapse = ""),
      test = if (!is.null(attr(d, "test"))) attr(d, "test") else if ("LRTPvalue" %in% colnames(mc)) "LRT" else if ("WaldPvalue" %in% colnames(mc)) "Wald" else NA,
      lrt_description = lrt_desc,
      full_model = sub("^LRT p-value: '(.*)' vs '(.*)'$", "\\1", lrt_desc),
      reduced_model = sub("^LRT p-value: '(.*)' vs '(.*)'$", "\\2", lrt_desc),
      reduced_matrix_cols = paste(colnames(attr(d, "reducedModelMatrix")), collapse = ";"),
      betaPrior = if (!is.null(attr(d, "betaPrior"))) attr(d, "betaPrior") else NA,
      dispersion_fitType = if (!is.null(dispersionFunction(d))) attr(dispersionFunction(d), "fitType") else NA,
      deseq2_version_in_object = as.character(metadata(d)$version),
      group_var = grp,
      group_levels = paste(levels(factor(cd[[grp]])), collapse = ";"),
      baseline = levels(factor(cd[[grp]]))[1],
      colData_cols = paste(colnames(cd), collapse = ";"),
      scaled_equals_scale_nCount = if (!is.null(d$scaled)) isTRUE(all.equal(as.numeric(d$scaled), as.numeric(scale(d$nCount_RNA, center = FALSE)))) else NA,
      cells_per_group = paste(names(table(cd[[grp]])), table(cd[[grp]]), sep = ":", collapse = ";"),
      cells_per_mouse = paste(names(table(cd$orig.ident)), table(cd$orig.ident), sep = ":", collapse = ";"),
      n_mice = length(unique(cd$orig.ident)),
      dup_mcols = paste(unique(colnames(mc)[duplicated(colnames(mc))]), collapse = ";"),
      assays = paste(assayNames(d), collapse = ";")),
      sf_provenance(d), filter_stats(counts(d)))
    mouse_cols <- intersect(c("orig.ident", "condition", "time", "sex", "age", "batch", "condition_time", "exposure", "type"), colnames(cd))
    mouse_tab <- unique(cd[, setdiff(mouse_cols, "type"), drop = FALSE])
    ranks <- list()
    if (isexp) {
      cd2 <- cd; cd2$condition_time[cd2$condition == "Naive"] <- "Naive"
      ranks$sbc <- design_rank(cd2, ~ sex + batch + condition)
      ranks$ssbc <- design_rank(cd2, ~ scaled + sex + batch + condition)
      ranks$sbct <- design_rank(cd2, ~ sex + batch + condition_time)
      ranks$sbe_cond <- design_rank(cd2, ~ sex + batch + condition)
      for (cc in c("Conditioned", "Unconditioned")) {
        s <- cd2[cd2$condition == cc, ]
        ranks[[paste0("time_", cc)]] <- design_rank(s, ~ sex + time)
        ranks[[paste0("time_batch_", cc)]] <- design_rank(s, ~ sex + batch + time)
      }
      tab <- as.data.frame(colData(d))
      rr <- extract_results(d, exp_contrasts)
    } else {
      ranks$stored_design <- design_rank(cd, design(d))
      if (all(c("sex", "batch") %in% colnames(cd))) ranks$sb_group <- design_rank(cd, as.formula(paste("~ sex + batch +", grp)))
      rr <- extract_results(d, list())
    }
    res$inventory <- inv; res$mouse_table <- mouse_tab; res$ranks <- ranks
    res$results <- rr; res$meta <- fit_meta(d, "stored")
    res$colData <- cd[, intersect(c("orig.ident", "condition", "time", "sex", "batch", "condition_time", "nCount_RNA", "type"), colnames(cd))]
    res$sizeFactors <- sizeFactors(d)
  } else {
    if (grepl("^exp_|^condtime_", analysis) && !isexp) stop("not an experience object")
    fit <- switch(analysis,
      exp_sheridan = { dd <- sheridan_wrap(d, ~ sex + batch + condition); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      exp_scaled   = { dd <- clean_from(d, ~ scaled + condition); wh(deseq_lrt(dd, ~ scaled, ncores)) },
      exp_sbc      = { dd <- clean_from(d, ~ sex + batch + condition); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      exp_ssbc     = { dd <- clean_from(d, ~ scaled + sex + batch + condition); wh(deseq_lrt(dd, ~ scaled + sex + batch, ncores)) },
      exp_scaled_N = { dd <- clean_from(d, ~ scaled + condition, relevel_naive = TRUE); wh(deseq_lrt(dd, ~ scaled, ncores)) },
      exp_sbc_N    = { dd <- clean_from(d, ~ sex + batch + condition, relevel_naive = TRUE); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      engram_sheridan = { e <- engram_cells_for(d); d$engram <- "no"; colData(d)[e, "engram"] <- "yes"
                          dd <- sheridan_wrap(d, ~ sex + batch + engram); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      engram_sbe   = { e <- engram_cells_for(d); d$engram <- "no"; colData(d)[e, "engram"] <- "yes"
                       dd <- clean_from(d, ~ sex + batch + engram); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      engram_sbce  = { e <- engram_cells_for(d); d$engram <- "no"; colData(d)[e, "engram"] <- "yes"
                       dd <- clean_from(d, ~ sex + batch + condition + engram); wh(deseq_lrt(dd, ~ sex + batch + condition, ncores)) },
      timeC_sheridan = { s <- d[, d$condition %in% "Conditioned"]; dd <- sheridan_wrap(s, ~ sex + time); wh(deseq_lrt(dd, ~ sex, ncores)) },
      timeU_sheridan = { s <- d[, d$condition %in% "Unconditioned"]; dd <- sheridan_wrap(s, ~ sex + time); wh(deseq_lrt(dd, ~ sex, ncores)) },
      timeC_clean  = { s <- d[, d$condition %in% "Conditioned"]; dd <- clean_from(s, ~ sex + time); wh(deseq_lrt(dd, ~ sex, ncores)) },
      timeU_clean  = { s <- d[, d$condition %in% "Unconditioned"]; dd <- clean_from(s, ~ sex + time); wh(deseq_lrt(dd, ~ sex, ncores)) },
      condtime_sheridan = { colData(d)[d$condition == "Naive", "condition_time"] <- "Naive"
                            dd <- sheridan_wrap(d, ~ sex + batch + condition_time); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      condtime_clean = { colData(d)[d$condition == "Naive", "condition_time"] <- "Naive"
                         dd <- clean_from(d, ~ sex + batch + condition_time); wh(deseq_lrt(dd, ~ sex + batch, ncores)) },
      # cell-type model on naive cells: ~ sex + batch + type is not full rank (the only female naive mouse, NC5, is
      # the only naive mouse in batch F), so we fit the closest estimable model, ~ batch + type vs ~ batch
      celltype_bt = { dd <- clean_from(d, ~ batch + type); wh(deseq_lrt(dd, ~ batch, ncores)) },
      stop("unknown analysis ", analysis))
    con <- if (grepl("^exp_", analysis)) exp_contrasts else
           if (grepl("^engram_", analysis)) list(EvN = c("engram", "yes", "no")) else
           if (grepl("^time", analysis)) list(T1v0 = c("time", "1hr", "0hr")) else list()
    res$results <- extract_results(fit, con)
    res$meta <- fit_meta(fit, analysis)
  }
  res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res$sessionInfo <- paste(capture.output(sessionInfo()), collapse = "\n")
  saveRDS(res, paste0(out_file, ".tmp"))
  file.rename(paste0(out_file, ".tmp"), out_file)
  message(format(Sys.time()), " done ", job, " in ", round(res$seconds), " s")
}
