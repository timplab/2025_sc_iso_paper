# DE model audit, GENOME/2025/281698 (27 Nov 2026)

**Case B.** All published DE results come from the saved fits: `~ scaled + condition` vs `~ scaled` (DESeq2 1.38.3 LRT;
scaled = nCount_RNA / RMS; C+ baseline; no sex or batch). This covers Supp Data 6, Fig 3D/S12/S15 and every number
quoted from them. The published lists match these fits exactly in 15/15 gene and 16/16 isoform cell types. Sheridan's
`~ sex + batch` calls never took effect, in either DESeq2 version. `DESeqDataSet()` on a fitted object relabels the
stored results as inputs, `DESeq()` appends the new fit under the same column names, and `results()` returns the
first, stored copy. Refitting the Methods model changes numbers and featured genes (1.46.0 shown. In 1.38.3,
the `scaled` model gives the same values and the literal model differs by at most 1 CA1 gene and 14 isoform rows;
table3b):

| | published | ~ sex + batch + cond | ~ scaled + sex + batch + cond |
|---|---|---|---|
| CA1 genes q < 0.05 | 597 | 288 | 324 |
| C+-up CA1 genes (also U+-up) | 250 (219, 87.6%) | 204 (198, 97.1%) | 255 (248, 97.3%) |
| Isoforms: rows / unique / genes | 1,972 / 1,266 / 1,217 | 4,869 / 2,429 / 2,297 | 1,526 / 995 / 969 |
| Median % isoforms, no gene change | 21.3 | 46.1 | 14.4 |
| Cell types with isoform changes (Olig merged) | 14 | 14 | 10 |
| Fos, Plp1 (CA1); Sgk1 (Olig) | C+U+, C+N; C+ | filtered, Cook's; C+U+ | ns, C+; C+U+ |
| Fos-201, Baiap2-204, Arpc2-201, Arpc2-203 | C+U+, C+U+, C+U+, U+ | ns, C+U+, ns, ns | ns, C+U+, ns, ns |

Syngap1-204 is not experience-variable under any model. The literal Methods model has no real depth normalization:
each isoform object carries a dummy all-1 transcript, so its size factors are 1. Gene size factors come from the 1-3
genes detected in every cell. That is why isoform calls inflate to 4,869. If we add sex and batch, keep `scaled`. All
designs are full rank (each mouse is one batch), and the baseline level does not change the labels.

**Time models (R2 #3, Fig S9).** These have the same problem. The five-level isoform list is identical to Supp Data 6.
The within-experience lists match the experience lists (median Jaccard 0.91-0.99). Clean refits share almost nothing
with them (median Jaccard 0.04 for genes, 0 for isoforms). Two copy-paste errors also make the Microglia and Olig
five-level gene rows copies of Inhibitory and Neuroblast. Regenerate before quoting.

**Supp Data 4** is `~ scaled + type` (8,908 isoforms, exact match). Sex and batch cannot both enter among naive mice.
`~ batch + type` gives 9,942 isoforms in both versions (Jaccard 0.88).

**Engram refit (A11).** Published Fig 4C/4D/S21B show experience-model q-values with engram fold changes (143 genes,
204 isoforms). Clean refits:

| | genes q<0.05 (up/down, abs LFC > 0.25) | isoforms |
|---|---|---|
| `~ sex + batch + engram` (1.38.3 / 1.46.0) | 36 (30/4) / 130 (80/27) | 592 (589/2) both |
| `~ sex + batch + condition + engram` (both) | 98 (60/17) | 587 (584/2) |

- **"Plp1 down in engram cells" does not survive.** Plp1 and Plp1-201 are Cook's outliers in every clean refit
  (maxCooks 8.7-15.9 vs cutoff 2.3-2.5), matching the three-droplet doublet result.
- **The other featured engram genes hold.** Junb, Egr1 and Fos-201, Egr1-201, Junb-201, Egr4-201, Arpc2-201 stay up.
  Fos and Egr4 lose their adjusted p to independent filtering in 3 of 4 refits (raw p 1.7e-151, 1.3e-6).
- **Panels 4C, 4D, 4E, S21B and S21C all change.** Fig 4E as coded forces Plp1 into the panel and orders by Fos, so it
  cannot be redrawn unchanged.
- **36 vs 130 is not Cook's filtering.** Cook's removes only 3 genes in either version. The cause is DESeq2 1.39.7's
  new independent-filtering fallback: 1.38.3 fell back to no filtering (36 = plain BH), while 1.46.0 filters 17,814
  low-count genes (130). The same rule change explains CA2 (72 vs 264 genes on the saved model).
- **Mixed, your call.** 1.46.0 is the corrected release but filters out Fos itself. The condition-adjusted model gives
  the same answer in both versions.

**Threshold (M24).** The code uses (mean + SD)/2 = 0.254 (z = 0.29), not z > 0.5. Both this and a true z > 0.5 (0.330)
select the same 150 cells, the CA1 cells with at least one Fos UMI (smallest nonzero 0.437). This reproduces
`engram_list.rds` exactly: C+ 85/408 (20.8%), U+ 58/333 (17.4%), N 7/102 (6.9%).

**Filters (M25).** Genes: detected in at least 10 cells overall, with no per-cell-type filter. Isoforms: no filter
(inferred; the object-building code was not found). Every annotated transcript is in the matrix, and 23,200-55,138
non-zero isoforms are tested per cell type. Table 6 is ready for Supp Data 5.
