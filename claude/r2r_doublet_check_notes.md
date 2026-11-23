# Notes for agents resuming the R3 doublet work (written 2026-11-23)

- /atium is mounted on zoro and winry, not duchess. Sheridan's R library (Seurat 5.3.0, DoubletFinder 2.0.6, DESeq2 1.38.3, R 4.2.2)
  is at winry:/home/sherbear/R/x86_64-pc-linux-gnu-library/4.2 (usable read-only via R_LIBS_USER).
- R2C2_combined_final.rds still contains NC1 (1,959 cells). Drop it: 16 mice = 17,800 cells, CA1 = 843 (408 C+, 333 U+, 102 N).
  Idents are already named ("Oligo. 1", "Gliaform" = OPC); Sheridan renames them by level order (see Rmd).
- The CA1 gene object genes_scaled_cond_nobatch.rds is a FITTED DESeqDataSet, design ~ scaled + condition, LRT vs ~ scaled,
  baseline Conditioned. A fresh refit (same size factors, fitType local, useT, minmu 1e-6, minReplicatesForReplace Inf) reproduces it
  exactly (597 LRT genes, 250 C+-up genes, Plp1 q 0.00147). The manuscript Methods say ~ sex + batch + experience with N baseline.
- Sheridan's before/after doublet comparison subset this fitted object and called results() without refitting: q-values unchanged by construction.
- Figure_4_plotting.Rmd calls DESeq() on that already-fitted object for the engram model. In DESeq2 1.38.3 this leaves duplicate
  mcols columns and results() returns the OLD experience-model LRT p-values (597 genes; identical p). A clean refit (mcols(d) <- NULL)
  gives 36 engram genes in 1.38.3 and 130 in 1.46.0; Plp1 is Cook's-filtered (maxCooks 11.7-13.7 vs cutoff 2.53).
- Engram threshold in code is (mean + sd)/2 of normalized Fos, not mean + 0.5 sd as the Fig 4A legend says.
- Key result: the C+ Plp1 increase in CA1 comes from 3 droplets (FC6 151 Plp1 UMIs, FC1 54 and 43) that also carry Mag/Mog/Mobp/Cldn11.
  Removing the top one: log2FC 1.06 -> 0.44, q 0.17; top three: -0.26, q 0.49. Dropping mouse FC6: 0.47, q 0.13.
- Sheridan's three DoubletFinder lists are not reproducible (pK not recorded) and agree poorly with a fresh per-sample run (Jaccard 0.08-0.17).
