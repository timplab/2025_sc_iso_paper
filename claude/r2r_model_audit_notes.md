# Notes for agents resuming the R2R DE model audit (written 2026-11-26)

Branch `r2r-model-audit` (from `r2r-doublets`). Code: `r2r_model_audit/`. Outputs: `/atium/Data/projects/sheridansinglecell/r2r_model_audit/`
(`cache/<DESeq2 version>/`, `tables/`, `source_data/`, `logs/`, `r2r_model_audit.html`, `sessionInfo.txt`, `code/` = copy of this folder).

- Hosts: /atium is on zoro and winry only. DESeq2 1.38.3 = Sheridan's R 4.2.2 library on winry
  (`R_LIBS_USER=/home/sherbear/R/x86_64-pc-linux-gnu-library/4.2`, `/usr/bin/Rscript`). DESeq2 1.46.0 = zoro conda env `r2r_doublets`.
  NFS caching: after editing code on /atium from zoro, winry may read the old file for a minute.
- Every saved per-cell-type DE object (genes_scaled_cond_nobatch, tx_scaled_cond_nobatch, tx_se_gencode_LRT, Olig_scaled_cond_nobatch)
  is a fitted DESeqDataSet: LRT `~ scaled + condition` vs `~ scaled`, baseline Conditioned, DESeq2 1.38.3, fitType local.
  Isoform objects have size factors all 1; gene objects have `type="ratio"` median-of-ratios size factors computed over the
  1-5 genes detected in every cell. tx_se_gencode_LRT objects hold all 149,384 GENCODE vM32 transcripts (no expression filter).
- The wrap bug: `DESeqDataSet(fitted, design=...)` relabels every existing rowData column as type "input", so `DESeq()` does not
  strip the old results; new results are appended with duplicate names and `results()` returns the first (stored) copy.
  Same code in 1.38.3 and 1.46.0. Affects DESeq2_analysis.Rmd (all experience lists), Figure_4 (engram genes), Figure_7 (engram
  isoforms), plotting_for_revision.Rmd (within-condition time) and winry:/home/sherbear/condition_time_for_supp.Rmd (five-level time).
  For the engram and time models the LFCs are new (new coefficient names) but the q-values are the stored experience LRT.
- condition_time_for_supp.Rmd also has copy-paste bugs: Microglia and Olig gene blocks lack `results()` (their rows in
  cond_time_sig_genes.tsv are copies of Inhibitory and Neuroblast); the Progenitor gene block uses `~ sex + batch + condition`.
- DESeq2 1.44.0 NEWS: LRT Cook's distances now use the full-model mean. This changes Cook's filtering for designs with replicated
  model-matrix rows (sex + batch + engram) but not for the saved `~ scaled + condition` fits: in CA2 and CA3 every maxCooks is far below the cutoff, and in CA1 maxCooks is NA.
- Independent filtering differs between versions. DESeq2 1.39.7 NEWS: when no filter threshold reaches the target, IF now falls
  back to the threshold giving 90%/80% of the maximum rejections; 1.38.3 falls back to no filtering. This explains engram
  `~ sex + batch + engram` genes 36 (1.38.3, = plain BH) vs 130 (1.46.0) and CA2 genes 72 vs 264 on the saved model.
- Engram threshold: coded (mean + SD)/2 = 0.2543 on log-normalized Fos (scale.factor = mean nCount_RNA = 4949.18); true z > 0.5 =
  0.3303; both select the same 150 cells, identical to engram_list.rds (Seurat 5.3.0). Seurat 5.5.1 fails on `slot=`.
- Google Docs (R2R, manuscript, supplement) need the user's approval card for read_doc; they were not read in this session.

## Status (2026-11-27)
All 368 fits cached for both versions; final render OK (`render.sh`, zoro). Verdict Case B. One-page summary and draft text are in
`r2r_model_audit/deliver/` and in the /atium output folder. Open items for Winston: choose the experience model (Option A or B),
choose the engram model (`~ sex + batch + engram` or `~ sex + batch + condition + engram`, and the DESeq2 version), confirm the Chawla 2013 PMID (23158763), approve or export the Google Docs so the
R2R text can be checked against the manuscript, regenerate the time lists before quoting them.
