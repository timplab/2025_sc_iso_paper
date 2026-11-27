# Draft text for the GENOME/2025/281698 revision (model audit, 2026-11-27)

Plain drafts for Winston to place. Numbers come from the audit tables in
`/atium/Data/projects/sheridansinglecell/r2r_model_audit/tables/`. Bracketed text marks a choice you need to make.
This is Case B, so two options are written out.
**Option A** keeps the published results and describes the model that was actually fitted.
**Option B** refits with sex and batch. If you choose B, we recommend `~ scaled + sex + batch + condition` over
`~ sex + batch + condition`. Without the `scaled` covariate, the saved size factors give no usable depth normalization
(see summary).

---

## 1. Methods, "Differential expression analysis"

### Option A (describe what was fitted; published numbers stand)

We tested for differential expression in each cell type separately, at the gene level and at the isoform level, with
DESeq2 v1.38.3 (Love et al. 2014, PMID 25516281). Each cell was one observation. We used the settings recommended for
single-cell counts (Van den Berge et al. 2018, PMID 29478411): a local dispersion trend (fitType = "local"), minmu =
1e-6, no replacement of outlier counts (minReplicatesForReplace = Inf), and a likelihood ratio test (LRT). To adjust for
sequencing depth we included each cell's total UMI count as a continuous covariate, scaled by its root mean square
across the cells of that type (R `scale()` without centering; "depth" below). Gene-level size factors were DESeq2
median-of-ratios estimates. Isoform-level size factors were set to 1, so depth enters the isoform models only through
the covariate. Sex and batch were not included in these models (see Reply to Reviewer 2, point 3).

For experience we compared the full model ~ depth + experience with the reduced model ~ depth. Experience is a
three-level factor (C+, U+, N). For cell type we fit the same pair of models to naive cells only, with cell type in
place of experience (~ depth + cell type vs ~ depth). The LRT compares full and reduced models, so its p-value does not
depend on which level is the reference. We used C+ as the reference. Using N instead left every feature's experience
label unchanged in CA1 [one low-count gene, Foxg1, changed in one of eight checks; table3e].

We called a feature experience-variable at a Benjamini-Hochberg adjusted LRT p-value below 0.05, after DESeq2's
default independent filtering on mean normalized count and default Cook's distance filtering. We assigned each
experience-variable feature to the experience(s) in which it was higher. For this we used the maximum-likelihood log2
fold changes from the same fit for the three pairwise comparisons (C+ vs U+, C+ vs N, U+ vs N). A feature was up in
an experience if its log2 fold change against at least one other experience was above 0.25. Features up in two
experiences relative to the third were labeled with both (C+U+, C+N, U+N). Features that passed the LRT but had no
pairwise log2 fold change above 0.25 were labeled "log2 < 0.25".

Before testing, we kept genes detected (at least one UMI) in at least 10 cells across the full dataset. We did not
apply an expression filter to isoforms before testing; DESeq2 excludes features with zero counts in every cell of the
cell type being tested. Supplementary Data 5 lists, for each cell type and level, the number of features in the
matrix, the number tested, the number with an adjusted p-value after independent and Cook's filtering, and the number
significant.

Engram analysis [same text for Option A and B; pick one model, see section 5]: we compared engram and non-engram CA1
cells with the LRT ~ sex + batch + engram vs ~ sex + batch [or ~ sex + batch + experience + engram vs ~ sex + batch +
experience, which accounts for the different proportions of engram cells across experiences]. We fit this with DESeq2
v[1.46.0] and the settings above, and called engram-variable features at adjusted p < 0.05 and |log2 fold change| >
0.25.

Time analysis [see section 4 for the formulas]: we compared ~ sex + time with ~ sex within C+ cells and within U+
cells, separately for each cell type (time = 10 or 60 min after exposure). Batch is not in this model because, within
an experience, time cannot be estimated once batch is included (the design is not full rank). [If you use the
five-level model instead: ~ sex + batch + condition_time vs ~ sex + batch, condition_time = C+ 10 min, C+ 60 min, U+ 10
min, U+ 60 min, N.] Note: the published time lists do not come from these models (see summary). They need to be
regenerated.

### Option B (refit with sex and batch; numbers change)

As in Option A, with these changes. For experience we compared ~ depth + sex + batch + experience with ~ depth + sex +
batch, fit with DESeq2 v[1.46.0]. Sex and batch were fixed effects in every model. Each mouse was processed in a
single batch (batches A to F; Table S4). Experience effects are therefore estimated from comparisons between mice of
different experiences within the same batch. All designs were full rank. For cell type, sex and batch are confounded
among the four naive mice: the only female naive mouse is also the only naive mouse in batch F. The cell-type model
therefore includes batch only: ~ depth + batch + cell type vs ~ depth + batch.

[If you prefer the literal current Methods model, ~ sex + batch + experience vs ~ sex + batch, say instead that depth
is handled only by size factors. Isoform size factors would then need to be recomputed first, because the saved
isoform objects have size factors fixed at 1.]

---

## 2. Fig 4A legend and Methods sentence for the engram threshold

**Fig 4A legend (replace the z > 0.5 clause):** Engram cells (red) are CA1 cells with detectable Fos: normalized Fos
above a threshold of (mean + SD)/2 across all CA1 cells, which selects every cell with at least one Fos UMI. This
calls 20.8% of C+ (85/408), 17.4% of U+ (58/333) and 6.9% of N (7/102) CA1 cells.

**Methods sentence:** We defined engram cells from log-normalized Fos expression in CA1 (Seurat LogNormalize, scale
factor equal to the mean CA1 UMI count, 4,949). We set the threshold at (mean + SD)/2 across all 843 CA1 cells (0.254,
which is 0.29 SD above the mean). We chose it because it called about 21% of C+ CA1 cells (85 of 408), close to the
roughly 25% expected from prior work (Chawla et al. 2013, PMID 23158763 [confirm this is the Chawla 2013 cited in
Results]). The smallest nonzero Fos value was 0.437, so this threshold selects exactly the 150 CA1 cells with at least
one Fos UMI. A threshold of mean + 0.5 SD (0.330) selects the same cells.

---

## 3. Reply to Reviewer 2, point 5 (isoform distribution and filters)

We thank the reviewer and agree that the filters should be stated explicitly. At the gene level we kept genes detected
in at least 10 cells across the full dataset before testing. At the isoform level we did not apply an expression
filter before testing: each cell-type matrix contained every annotated transcript, and DESeq2 tested the transcripts
with at least one count in that cell type. DESeq2 then removed low-count features by independent filtering and
outliers by Cook's distance before multiple-testing correction. Supplementary Data 5 now reports, for each cell type,
the number of genes and isoforms in the matrix, the number tested, the number retained after filtering (with an
adjusted p-value), and the number significant. For example, in CA1 we tested 19,552 genes; 4,010 had an adjusted
p-value and 597 were significant. In CA1 we also tested 39,276 isoforms; 5,888 had an adjusted p-value and 284 were
significant. Across cell types we tested 23,200 to 55,138 isoforms. [Option B: replace the CA1 numbers with the refit
values from table6/table3.]

---

## 4. Reviewer 2, point 3: time-model formulas for the brackets

[within each experience, separately for C+ and U+ cells of each cell type: ~ sex + time vs ~ sex, DESeq2 LRT; batch
is not included because time is not estimable within an experience once batch is in the model]

Alternative, if Fig S9 uses the five-level model: [~ sex + batch + condition_time vs ~ sex + batch, with
condition_time = C+ 10 min, C+ 60 min, U+ 10 min, U+ 60 min, N]

Caution: the published time lists (time_sig_genes.tsv, time_sig_isos.tsv, cond_time_sig_genes.tsv,
cond_time_sig_isoforms.tsv) carry the stored experience q-values, not these models. The reply should quote numbers
from the clean refits (table2, rows "Within-condition time model" and "Five-level condition_time model", models
clean_*).

---

## 5. Replacement Results sentences for the engram paragraphs

The published Fig 4C/4D and Fig S21B q-values are the experience-model q-values, so these paragraphs change under any
valid engram model. Two versions follow. Numbers are from DESeq2 1.46.0.

### 5a. Model ~ sex + batch + engram (the model as coded, refit cleanly)

Fig 4C/4E paragraph: Comparing engram with non-engram CA1 cells, we found 130 differentially expressed genes (adjusted
p < 0.05). Of these, 107 had |log2 fold change| > 0.25 (80 higher and 27 lower in engram cells). As expected,
immediate early genes were higher in engram cells, including Junb and Egr1. Fos defines the engram cells, so it is
higher by construction; DESeq2's independent filtering removed it, and Egr4, from the adjusted test because of their
low mean counts (unadjusted LRT p = 1.7e-151 for Fos and 1.3e-6 for Egr4). Plp1 was not differential: a few outlying cells drove its
apparent decrease, and DESeq2's Cook's distance filter removed it (see also the doublet analysis). [Fig 4E sentence:
rebuild the co-expression panel from the new list. As coded, the panel is ordered by correlation with Fos and always
includes Plp1, so it cannot be drawn unchanged.]

Fig S21 paragraph: At the isoform level, 592 isoforms differed between engram and non-engram CA1 cells (adjusted p <
0.05). Of these, 591 had |log2 fold change| > 0.25, nearly all higher in engram cells (589 up, 2 down). They included
Fos-201, Egr1-201, Junb-201 and Egr4-201. Plp1-201 was removed by Cook's distance filtering and is not differential.

### 5b. Model ~ sex + batch + experience + engram (adjusts for the uneven share of engram cells across experiences)

Fig 4C/4E paragraph: After accounting for experience, 98 genes differed between engram and non-engram CA1 cells
(adjusted p < 0.05). Of these, 77 had |log2 fold change| > 0.25 (60 higher and 17 lower in engram cells), including
Junb and Egr1. As in 5a, independent filtering removed Fos and Egr4 from the adjusted test (unadjusted LRT p = 1.8e-149
and 4.7e-6). Plp1 was removed by Cook's distance filtering and is not differential. [Same Fig 4E caveat as 5a.]

Fig S21 paragraph: At the isoform level, 587 isoforms differed. Of these, 586 had |log2 fold change| > 0.25 (584
higher and 2 lower in engram cells), including Fos-201, Egr1-201, Junb-201 and Egr4-201. Plp1-201 was removed by
Cook's distance filtering.

The counts for 5b, and for the isoform model in 5a, are the same in DESeq2 1.38.3 and 1.46.0. Only the 5a gene count
depends on the version (36 genes in 1.38.3; see summary).
