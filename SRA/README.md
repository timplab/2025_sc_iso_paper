# SRA submission: PacBio Kinnex data → PRJNA1368682

Add Sheridan's PacBio Kinnex single-cell reads to the BioProject that already
holds the ONT R2C2 data for "Single-cell long-read sequencing of the
experience-induced transcriptome".

## What already exists

| | |
|---|---|
| BioProject | **PRJNA1368682** (SRA study SRP648294), public since 2025-11-24 |
| Original SRA submission | SUB15610104, submitted under Winston Timp's NCBI account (wtimp1) |
| BioSamples | 17, `2025_SC_{FC1-6,UC1-6,NC1-5}` = SAMN53367189–SAMN53367205 |
| Runs | 17 ONT PromethION fastq.gz runs, library_IDs `SC_<sample>_t<min>` |

PRJNA1327322 has the same title but no data attached. It is an earlier
duplicate and should not be used.

## Data to submit

`/atium/Data/projects/sheridansinglecell/pacbio/261006_sherbear_pacbio_kinnex_bams/`:
16 BAMs, 144 GB.

- These are the instrument's unaligned HiFi read BAMs (Revio, ICS 13.5, on-instrument
  CCS + jasmine + lima Kinnex barcode demux, barcodes bcM0001–4). Each read is
  an unsegmented Kinnex 16-mer array (FC1: 2.5M reads → 39.9M S-reads, mean
  concatenation 15.8).
- 4 SMRT cells (pools 1_rpt, 2, 3, 4) × 4 barcodes = 16 samples. Each sample is on one
  SMRT cell, so the submission has one BAM → one SRA run → one experiment per sample.
- The `@RG SM` tag in every BAM matches the sample prefix of its filename (checked).
- There is no PacBio data for NC1.
- **NC2 (60 MB) and NC4 (30 MB, 6,257 reads) are far smaller than the others** (1–15 GB each).

We submit these raw HiFi BAMs, not the skera-segmented or Iso-Seq-processed BAMs in
`pacbio/<sample>/`. SRA prefers the least-processed form, and anyone can redo
segmentation with `skera split`. We should note this in the design description,
which the draft does.

## Decisions to make before submitting

1. **Same BioSamples?** The draft assumes the Kinnex libraries were made from the
   same 10X cDNA and animals as the ONT runs. If so, the new runs attach to the
   existing SAMN accessions and no new BioSamples are needed. If they came from new
   animals or tissue, we need new BioSamples, and the metadata sheet changes.
2. **NC2 / NC4**: submit them anyway (recommended if they are used anywhere in the
   paper), or drop them.
3. **pool_1**: FC1, FC5, NC5 and UC4 come from a repeat run (`pool_1_rpt`). Was there
   an original pool_1 SMRT cell? If it has usable data, it should be a second run
   under the same experiments.
4. **Release date**: immediately (the project is already public), or hold until
   publication.

## Steps

1. **Pre-flight check** (quickcheck, check SM tag against library_ID, md5):
   ```
   ./check_bams.sh /atium/Data/projects/sheridansinglecell/pacbio/261006_sherbear_pacbio_kinnex_bams \
       /home/timp/miniforge3/envs/240811_bcmb/bin/samtools
   ```
   This writes `md5sums.txt`. md5 over 144 GB on /atium takes a while, so run it in tmux.

2. **Review `pacbio_sra_metadata.tsv`**: a draft SRA metadata sheet with one row per BAM.
   Fields match the existing ONT experiments (RNA-Seq / TRANSCRIPTOMIC SINGLE CELL /
   PolyA / single) except `platform=PACBIO_SMRT`, `instrument_model=Revio`,
   `filetype=bam`. Library IDs are `SC_<sample>_t<min>_PB`, using the existing
   timepoints, so they can't collide with the ONT library IDs. Check the design
   description wording, then paste the sheet into the portal's SRA metadata
   template (the portal checks the column names and allowed values).

3. **Upload** (needs the Aspera key and preload folder from the portal: *My Preload
   Folder → Aspera command line*). `ascp` is already in `~/bin`.
   ```
   ./upload.sh ~/aspera.openssh uploads/<account_folder>
   ```
   Files go to `<preload>/pacbio_kinnex/`. Expect a few hours for 144 GB. Wait
   about 10 minutes after the transfer finishes before selecting the folder in
   the portal.

4. **Create the SRA submission** at https://submit.ncbi.nlm.nih.gov/subs/sra/,
   logged in as the account that owns PRJNA1368682 (wtimp1), or one added to that
   submission group.
   - BioProject: *Yes, existing* → `PRJNA1368682`
   - BioSample: *Yes, existing* → the accessions are in the metadata sheet
   - SRA metadata: upload `pacbio_sra_metadata.tsv`
   - Files: *FTP or Aspera Command Line preload* → select `pacbio_kinnex`
   - Release date as decided above → Submit

5. **After SRA processes the submission** (usually a day or two):
   - Check that 16 new SRR runs show up on SRP648294 and that their spot counts
     match the BAM read counts.
   - Save the SUB/SRX/SRR accessions to `pacbio_sra_accessions.tsv` here.
   - Update the paper's data availability statement and the top-level README.

## Related cleanup on the existing records

- **BioProject description** mentions only ONT/R2C2. Add a sentence on the
  PacBio Kinnex (Revio) data. Edit it in the portal or email bioprojecthelp@ncbi.nlm.nih.gov.
- **BioSample titles are broken.** All 17 have an Excel formula as their title
  (`=RIGHT(A5,LEN(A5) - (FIND(CHAR(160),...`) instead of text. Fix it in the same
  pass: email biosamplehelp@ncbi.nlm.nih.gov with a SAMN → title table, or use
  the portal's BioSample update.

## Files here

| File | Purpose |
|---|---|
| `README.md` | This plan |
| `pacbio_sra_metadata.tsv` | Draft SRA metadata, one row per BAM |
| `check_bams.sh` | Pre-flight integrity, SM-tag, and md5 checks |
| `upload.sh` | Aspera upload of the files listed in the metadata sheet |
