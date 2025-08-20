# miRNA-Seq - basic analysis pipeline (Ion Torrent) — README

A simple, end-to-end **R** pipeline for miRNA differential expression analysis from Ion Torrent exports.  
It builds a unified count matrix, performs QC, normalization, descriptive DE (or DESeq2 if replicates exist), generates standard plots, and exports an integrated log2FC heatmap across contrasts.

To note: was developed for an A549 LC cell lines sequencing experiment, disregard hard-coding related to this research question. 
---

## Features

- Auto-ingests multiple CSVs and constructs a **counts matrix** (miRNA × samples).
- **QC**: library sizes, sample–sample correlation heatmap, PCA (VST).
- **Filtering** (low counts), **normalization** (edgeR TMM; DESeq2 VST/rlog).
- **Contrasts**:
  - Works **descriptively** when there are no replicates (log2FC only).
  - Uses **DESeq2** statistics automatically if each group has ≥2 replicates.
- **Plots**: volcano (safe fallback when no p-values), top5 bar & lollipop, heatmap (top variable miRNA), global log2FC heatmap across all contrasts.
- **Threshold summaries** at common log2FC cutoffs (≈1.5×, 2×, 2.5×).

---

## Requirements

- R ≥ 4.1
- Packages (auto-install on first run):  
  `tidyverse, data.table, DESeq2, edgeR, pheatmap, RColorBrewer, EnhancedVolcano, ggrepel, vsn`

> If your environment blocks installs, pre-install these packages or set a CRAN mirror accessible to your machine.

---

## Quick start

1) **Place your CSV files** in a folder (default expected path in the script):
```
analiza_zip/analiza/
  IonXpressRNA_011_rawlib.basecaller.bam [CARBOPLATIN] (single) Small RNA sample grouped on mature.csv
  IonXpressRNA_010_rawlib.basecaller.bam [VINORELBIN] (single) Small RNA sample grouped.csv
  ...
```

2) **Edit CONFIG** in `pipeline_miRNA_A549.R` if needed:
- `INPUT_DIR` (where your CSVs are)
- `OUT_DIR` (where results go)
- `label_to_group()` mapping (how file labels map to experimental groups)
- `CONTRASTS` (pairs to compare)

3) **Run**:
```bash
Rscript pipeline_miRNA_A549.R
```

---

## Input expectations

- **“grouped on mature” CSVs** containing columns:
  - `Name` (miRNA symbol; e.g., `hsa-miR-21-5p` / `let-7a-1`…)
  - `Total` (preferred) or `Expression values` (fallback) as counts
- The **sample label** is parsed from the filename text in square brackets `[ ... ]` by default (adjust `extract_label()` if needed).

---

## Outputs (folder structure)

```
results_A549/
  sample_metadata_auto.csv         # inferred sample→group mapping
  counts_filtered.csv              # counts after low-expression filtering

  qc/
    library_sizes.png
    correlation_heatmap.png
    PCA_vst.png

  normalized/
    normalized_counts.csv
    vst_matrix.csv                 # if available

  de/
    <contrast>_<grp1>_vs_<grp2>.csv   # log2FC (+ pval/padj if replicates)
    log2FC_all_contrasts.csv          # integrated matrix across contrasts
    counts_log2FC_ge_0.585.csv        # ~1.5× threshold summary
    counts_log2FC_ge_1.000.csv        # 2× threshold summary
    counts_log2FC_ge_1.322.csv        # ~2.5× threshold summary

  plots/
    volcano_<contrast>.png
    bar_top5_<contrast>.png
    lollipop_top5_<contrast>.png
    global_heatmap_top30.png
```

---

## Configuration notes

- **Group mapping**  
  Adjust `label_to_group()` to map your bracket labels to canonical group names.

- **Contrasts**  
  Define as `c("GroupA", "GroupB")`. The script also provides a **collapsed “VINO+CARBO” pseudo-group** by summing VINO and CARBO counts if you lack replicates—useful for descriptive comparisons. For proper inference, supply true replicates per condition.

- **No replicates?**  
  The pipeline will:
  - Skip DESeq2 model fitting.
  - Report **descriptive log2FC** only (no valid p-values/FDR).
  - Still produce volcano-style visuals with a **pseudo score**, clearly for visualization only.

- **With replicates (≥2/arm)**  
  DESeq2 runs automatically; results include **p-values and FDR**. You can extend to a full **factorial design** if you edit the `design` and result extraction.

---

## Interpretation caveats

- **Without biological replicates**, statistical inference is not valid. Treat outputs as **descriptive**.
- For publication-grade analyses: include **≥3 replicates** per arm and consider a **factorial DESeq2 model**.

---

## Customization tips

- **miRNA identifiers**: If you use MIMAT accessions instead of `Name`, adjust the column used in `read_one()`.
- **Filtering**: tune `keep <- rowSums(cpm > 1) >= …` to your sample size.
- **Top N plots**: set `TOP_N_VAR` (default 30) and `VOLCANO_TOP_N` (default 5).
- **Threshold summaries**: edit `LOG2FC_THRESHES` to your preferred cutoffs.

---

## Reproducibility

- Record R session info in your run logs:
  ```r
  sessionInfo()
  ```
- Pin package versions in a **`renv`** or **Docker** image for stable builds.

---

## Troubleshooting

- **Empty/garbled plots**: ensure the CSVs contain `Name` and `Total`/`Expression values`.
- **No PCA**: VST may fail if inputs are too small; QC still outputs library sizes and correlations.
- **No DESeq2 stats**: check `table(group)` in `sample_metadata_auto.csv`. You need ≥2 replicates per group.

---

## License & citation

- License: MIT.
- If this pipeline supports a publication, cite **DESeq2** and **edgeR** appropriately.

---

## Contact

- Maintainer: *Stefan Strilciuc*  
- Issues: open a GitHub issue with **sessionInfo()**, example filenames, and a short log excerpt.

---

**Happy analyzing!**
