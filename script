#!/usr/bin/env Rscript

# =============================
#  A549 miRNA-Seq pipeline (Ion Torrent / CLC exports)
#  Author: <you>
#  Version: 1.0
#  Requires: R >= 4.1
# =============================

suppressPackageStartupMessages({
  pkgs <- c(
    "tidyverse","data.table","DESeq2","edgeR","pheatmap",
    "RColorBrewer","EnhancedVolcano","ggrepel","vsn"
  )
  miss <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(miss) > 0) {
    message("Installing missing packages: ", paste(miss, collapse=", "))
    install.packages(miss, repos = "https://cloud.r-project.org")
  }
  lapply(pkgs, require, character.only = TRUE)
})

# ============== CONFIG =================
# Setează aici folderele și mapping-ul de contraste/grupuri

INPUT_DIR  <- "analiza_zip/analiza"   # directorul cu CSV-urile “grouped … .csv”
OUT_DIR    <- "results_A549"          # unde scriem rezultatele

# Heuristică pentru denumiri (editează dacă folosești alte label-uri în numele fișierelor)
# Eticheta probei este textul din paranteze [ ... ] din numele fișierului.
# Exemplu: "IonXpressRNA_011_rawlib.basecaller.bam [CARBOPLATIN] (single) Small RNA sample grouped on mature.csv"
#          -> label = "CARBOPLATIN"
label_to_group <- function(lbl) {
  lbl <- tolower(lbl)
  lbl <- gsub("\\s+", " ", lbl)
  case_when(
    grepl("negativ ctr inhibitor", lbl) ~ "CTRL_INH",
    grepl("negativ ctr mimic", lbl)     ~ "CTRL_MIM",
    grepl("^21 inh(\\b|\\s)", lbl)      ~ "MIR21_INH",
    grepl("181 mimic$", lbl)            ~ "MIR181_MIM",
    grepl("^vinorelb", lbl)             ~ "VINO",
    grepl("^carboplat", lbl)            ~ "CARBO",
    grepl("^21-v-c$", lbl)              ~ "MIR21_INH_V_C",
    grepl("^181-v-c$", lbl)             ~ "MIR181_MIM_V_C",
    grepl("^v si c$", lbl)              ~ "MIR21_INH_VINO_CARBO",
    grepl("^21-181-v-c$", lbl)          ~ "MIR21_INH_MIR181_MIM_V_C",
    TRUE                                ~ toupper(gsub("[^A-Za-z0-9]+","_",lbl))
  )
}

# Contraste (editează după cum ai designul)
CONTRASTS <- list(
  miR21_inh_vs_ctrl_inh             = c("MIR21_INH", "CTRL_INH"),
  miR181_mimic_vs_ctrl_mim          = c("MIR181_MIM", "CTRL_MIM"),
  miR21_inh__vino__carbo_vs_vinoCarbo = c("MIR21_INH_VINO_CARBO", "VINO+CARBO"),
  miR181_mimic__vino__carbo_vs_vinoCarbo = c("MIR181_MIM_V_C", "VINO+CARBO"),
  double_inh_mimic__vino__carbo_vs_vinoCarbo = c("MIR21_INH_MIR181_MIM_V_C", "VINO+CARBO")
)
# Observații:
# - Pentru comparațiile “... vs VINO+CARBO”, vom agrega Vino și Carbo simple ca un pseudo-grup “VINO+CARBO”
#   (util când nu ai replici). Dacă ai replici reale pe grupuri distincte, recomand design factorial (vezi mai jos).

# Praguri pentru vizualizări
TOP_N_VAR       <- 30
VOLCANO_TOP_N   <- 5            # etichetare doar top 5 up & 5 down
LOG2FC_THRESHES <- c(0.585, 1, log2(2.5))  # ~1.5x, 2x, 2.5x

# ============ UTILITARE I/O =============
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DIR_QC     <- file.path(OUT_DIR, "qc");         dir.create(DIR_QC, showWarnings = FALSE)
DIR_NORM   <- file.path(OUT_DIR, "normalized"); dir.create(DIR_NORM, showWarnings = FALSE)
DIR_DE     <- file.path(OUT_DIR, "de");         dir.create(DIR_DE, showWarnings = FALSE)
DIR_PLOTS  <- file.path(OUT_DIR, "plots");      dir.create(DIR_PLOTS, showWarnings = FALSE)

# ============ 1) CITIRE & MATRICE CONTURI ============
message(">>> Building count matrix from CSV files in: ", INPUT_DIR)
csvs <- list.files(INPUT_DIR, pattern = "grouped.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
stopifnot(length(csvs) > 0)

read_one <- function(f) {
  DT <- suppressWarnings(fread(f))
  # căutăm coloana de numărări: preferăm "Total", dacă nu există, folosim "Expression values"
  cols <- colnames(DT)
  count_col <- if ("Total" %in% cols) "Total" else if ("Expression values" %in% cols) "Expression values" else NA
  stopifnot(!is.na(count_col))
  stopifnot("Name" %in% cols)
  DT <- DT[, .(Name = as.character(Name), counts = as.integer(get(count_col))),]
  # agregare dacă există duplicate de "Name"
  DT <- DT[, .(counts = sum(counts, na.rm = TRUE)), by = Name]
  DT[]
}

extract_label <- function(path) {
  nm <- basename(path)
  m <- regexpr("\\[(.*?)\\]", nm)
  if (m[1] > 0) {
    lbl <- regmatches(nm, m)
    lbl <- gsub("^\\[|\\]$","", lbl)
  } else {
    lbl <- tools::file_path_sans_ext(nm)
  }
  lbl
}

# Colectăm
mat_list <- list()
labels   <- c()
for (f in csvs) {
  DT <- read_one(f)
  lbl <- extract_label(f)
  mat_list[[lbl]] <- DT
  labels <- c(labels, lbl)
}

# Union all names
all_names <- unique(unlist(lapply(mat_list, function(x) x$Name)))
mat <- matrix(0L, nrow = length(all_names), ncol = length(mat_list),
              dimnames = list(all_names, names(mat_list)))
for (lbl in names(mat_list)) {
  idx <- match(mat_list[[lbl]]$Name, all_names)
  mat[idx, lbl] <- mat_list[[lbl]]$counts
}
counts <- as.data.frame(mat, check.names = FALSE)

# ============ 2) METADATA (auto + manual) ============
sample_table <- tibble::tibble(
  sample  = colnames(counts),
  group   = sapply(colnames(counts), label_to_group)
)

# Construim un pseudo-grup “VINO+CARBO” (agregare) doar pentru comparațiile fără replici:
sample_table$group_collapsed <- dplyr::case_when(
  sample_table$group %in% c("VINO","CARBO") ~ "VINO+CARBO",
  TRUE                                      ~ sample_table$group
)

# Export pentru verificare/editare manuală
write_csv(sample_table, file.path(OUT_DIR, "sample_metadata_auto.csv"))

message("Sample table:")
print(sample_table)

# ============ 3) QC de bază ============
lib_sizes <- colSums(counts)
png(file.path(DIR_QC, "library_sizes.png"), width=1600, height=900, res=150)
barplot(sort(lib_sizes, decreasing = TRUE), las=2, main="Library sizes (raw counts)", ylab="reads")
dev.off()

# CPM + corelații
dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge)
cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
cor_mat <- cor(cpm, method = "pearson")
png(file.path(DIR_QC, "correlation_heatmap.png"), width=1200, height=1000, res=150)
pheatmap::pheatmap(cor_mat, main = "Sample–sample correlation (log-CPM)")
dev.off()

# PCA pe VST (dacă posibil)
vst_mat <- NULL
try({
  dds0 <- DESeqDataSetFromMatrix(countData = counts,
                                 colData   = data.frame(row.names=colnames(counts), group=sample_table$group),
                                 design    = ~ 1)
  vst_mat <- assay(vst(dds0, blind = TRUE))
  pca <- prcomp(t(vst_mat), scale. = FALSE)
  pc <- as.data.frame(pca$x[,1:2])
  pc$sample <- rownames(pc)
  png(file.path(DIR_QC, "PCA_vst.png"), width=1200, height=1000, res=150)
  ggplot(pc, aes(PC1, PC2, label = sample)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggtitle("PCA (VST)") +
    theme_minimal() -> g
  print(g)
  dev.off()
}, silent = TRUE)

# ============ 4) Filtrare expresie scăzută ============
keep <- rowSums(edgeR::cpm(dge) > 1) >= max(1, floor(ncol(counts) * 0.25))
counts_f <- counts[keep, , drop=FALSE]
write_csv(as.data.frame(counts_f), file.path(OUT_DIR, "counts_filtered.csv"))

# ============ 5) Normalizare & export ============
dds <- DESeqDataSetFromMatrix(countData = round(counts_f),
                              colData   = data.frame(row.names=colnames(counts_f), group=sample_table$group),
                              design    = ~ group)

# detectăm replici
has_reps <- all(table(colData(dds)$group) >= 2)

norm_counts <- NULL
vst_norm <- NULL
rlog_norm <- NULL

if (has_reps) {
  message("Replicates detected. Running DESeq2...")
  dds <- DESeq(dds, parallel = FALSE)
  vst_norm <- assay(vst(dds, blind = FALSE))
  rlog_norm <- assay(rlog(dds, blind = FALSE))
  norm_counts <- counts(dds, normalized = TRUE)
} else {
  message("No replicates per group. Skipping DESeq2 fit. Using edgeR TMM norm counts + VST blind.")
  dge2 <- edgeR::DGEList(counts = counts_f)
  dge2 <- edgeR::calcNormFactors(dge2)
  norm_counts <- edgeR::cpm(dge2, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0.1)
  # VST blind pe design ~1 (fără DE)
  dds0 <- DESeqDataSetFromMatrix(countData = round(counts_f),
                                 colData   = data.frame(row.names=colnames(counts_f), group="one"),
                                 design    = ~ 1)
  vst_norm <- assay(vst(dds0, blind = TRUE))
  rlog_norm <- assay(rlog(dds0, blind = TRUE))
}

write_csv(as.data.frame(norm_counts), file.path(DIR_NORM, "normalized_counts.csv"))
if (!is.null(vst_norm)) write_csv(as.data.frame(vst_norm), file.path(DIR_NORM, "vst_matrix.csv"))

# ============ 6) Funcții DE & plotting ============
safe_volcano <- function(df, title, out_png) {
  # df: data.frame cu col: miRNA, log2FC, pval, padj (pval/padj pot lipsi)
  if (!"padj" %in% names(df)) df$padj <- NA_real_
  if (!"pval" %in% names(df)) df$pval <- NA_real_
  df$label <- ""
  up   <- df %>% arrange(desc(log2FC)) %>% head(VOLCANO_TOP_N)
  down <- df %>% arrange(log2FC)       %>% head(VOLCANO_TOP_N)
  df$label[match(up$miRNA, df$miRNA)]   <- up$miRNA
  df$label[match(down$miRNA, df$miRNA)] <- down$miRNA

  png(out_png, width=1200, height=1000, res=150)
  if (requireNamespace("EnhancedVolcano", quietly=TRUE) && all(!is.na(df$pval))) {
    EnhancedVolcano::EnhancedVolcano(df,
      lab = df$label, x = "log2FC", y = "pval",
      title = title, pCutoff = 0.05, FCcutoff = 1,
      xlab = "log2FC", ylab = "-log10(p)",
      drawConnectors = TRUE, max.overlaps = Inf, legendPosition = "right")
  } else {
    # fallback ggplot (pseudop)
    df$neglog10p <- -log10(ifelse(is.na(df$pval), 1/(abs(df$log2FC)+1), df$pval))
    ggplot(df, aes(log2FC, neglog10p)) +
      geom_point(alpha=.6) +
      geom_vline(xintercept = 0, linetype="dashed") +
      ggrepel::geom_text_repel(aes(label = label), min.segment.length = 0, box.padding = .3, max.overlaps = Inf) +
      labs(title = title, x = "log2FC", y = "-log10(pseudo-p)") +
      theme_minimal()
  }
  dev.off()
}

bar_lollipop_top5 <- function(df, title, out_bar, out_lolli) {
  top_up   <- df %>% arrange(desc(log2FC)) %>% head(VOLCANO_TOP_N)
  top_down <- df %>% arrange(log2FC)       %>% head(VOLCANO_TOP_N)
  top <- bind_rows(top_up, top_down)
  # bar
  png(out_bar, width=1200, height=1000, res=150)
  ggplot(top, aes(x = reorder(miRNA, log2FC), y = log2FC, fill = log2FC > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("FALSE"="#377eb8","TRUE"="#e41a1c"), guide="none") +
    labs(title = paste("Top5 up/down –", title), x = "", y = "log2FC") +
    theme_minimal()
  dev.off()
  # lollipop
  png(out_lolli, width=1200, height=1000, res=150)
  ggplot(top, aes(x = log2FC, y = reorder(miRNA, log2FC), color = log2FC > 0)) +
    geom_segment(aes(x=0, xend=log2FC, y=reorder(miRNA, log2FC), yend=reorder(miRNA, log2FC))) +
    geom_point(size=3) +
    scale_color_manual(values = c("FALSE"="#377eb8","TRUE"="#e41a1c"), guide="none") +
    labs(title = paste("Lollipop –", title), x = "log2FC", y = "") +
    theme_minimal()
  dev.off()
}

heatmap_topn <- function(mat, title, out_png, n = TOP_N_VAR) {
  var_genes <- order(matrixStats::rowVars(as.matrix(mat)), decreasing = TRUE)[seq_len(min(n, nrow(mat)))]
  sel <- mat[var_genes, , drop = FALSE]
  # z-score pe rând
  sel_z <- t(scale(t(sel)))
  png(out_png, width=1400, height=1200, res=150)
  pheatmap::pheatmap(sel_z, main = title, show_rownames = TRUE, show_colnames = TRUE,
                     color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(255))
  dev.off()
}

# ============ 7) DE pe contraste ============
# pentru comparațiile cu “VINO+CARBO” agregăm VINO și CARBO dacă există ca probe separate
collapse_vc <- function(cnts, sTable) {
  vc_cols <- which(sTable$group %in% c("VINO","CARBO"))
  if (length(vc_cols) > 0) {
    new_col <- rowSums(cnts[, vc_cols, drop=FALSE])
    colnames_add <- paste0("VINO+CARBO_", length(vc_cols), "samp")
    cnts2 <- cbind(cnts, setNames(as.data.frame(new_col), colnames_add))
    sTab2 <- rbind(sTable, data.frame(sample = colnames_add, group = "VINO+CARBO", group_collapsed = "VINO+CARBO"))
    rownames(sTab2) <- sTab2$sample
    list(counts = cnts2, sample_table = sTab2)
  } else {
    list(counts = cnts, sample_table = sTable)
  }
}

coll <- collapse_vc(counts_f, sample_table)
counts_for_de <- coll$counts
sTable_de     <- coll$sample_table
stopifnot(all(colnames(counts_for_de) == sTable_de$sample))

# Helper: DE cu DESeq2 dacă are replici; altfel log2FC descriptiv
run_contrast <- function(cnts, sTab, grp1, grp2, out_prefix) {
  # selectăm coloanele
  cols1 <- which(sTab$group_collapsed == grp1 | sTab$group == grp1)
  cols2 <- which(sTab$group_collapsed == grp2 | sTab$group == grp2)
  if (length(cols1)==0 || length(cols2)==0) return(NULL)

  cnt1 <- cnts[, cols1, drop=FALSE]
  cnt2 <- cnts[, cols2, drop=FALSE]
  # log2FC descriptiv
  m1 <- if (ncol(cnt1)>1) rowMeans(cnt1) else as.numeric(cnt1)
  m2 <- if (ncol(cnt2)>1) rowMeans(cnt2) else as.numeric(cnt2)
  log2fc <- log2((m1+1)/(m2+1))

  res_df <- data.frame(miRNA = rownames(cnts), log2FC = log2fc, stringsAsFactors = FALSE)

  # Dacă am replici pe ambele brațe, încercăm DESeq2
  has_rep <- (ncol(cnt1) >= 2 && ncol(cnt2) >= 2)
  if (has_rep) {
    sub_counts <- cbind(cnt1, cnt2)
    grp <- c(rep(grp1, ncol(cnt1)), rep(grp2, ncol(cnt2)))
    dds_sub <- DESeqDataSetFromMatrix(countData = round(sub_counts),
                                      colData = data.frame(row.names=colnames(sub_counts), group = grp),
                                      design = ~ group)
    dds_sub <- DESeq(dds_sub)
    res <- results(dds_sub, contrast = c("group", grp1, grp2))
    res_df$pval <- res$pvalue
    res_df$padj <- res$padj
  }

  # scriem
  out_csv <- file.path(DIR_DE, paste0(out_prefix, "_", grp1, "_vs_", grp2, ".csv"))
  write_csv(res_df, out_csv)

  # ploturi
  safe_volcano(res_df, paste(grp1, "vs", grp2), file.path(DIR_PLOTS, paste0("volcano_", out_prefix, "_", grp1, "_vs_", grp2, ".png")))
  bar_lollipop_top5(res_df, paste(grp1, "vs", grp2),
                    file.path(DIR_PLOTS, paste0("bar_top5_", out_prefix, "_", grp1, "_vs_", grp2, ".png")),
                    file.path(DIR_PLOTS, paste0("lollipop_top5_", out_prefix, "_", grp1, "_vs_", grp2, ".png")))
  invisible(res_df)
}

all_res <- list()
for (nm in names(CONTRASTS)) {
  grp1 <- CONTRASTS[[nm]][1]
  grp2 <- CONTRASTS[[nm]][2]
  cat("Running contrast:", nm, " => ", grp1, "vs", grp2, "\n")
  r <- run_contrast(counts_for_de, sTable_de, grp1, grp2, out_prefix = nm)
  all_res[[nm]] <- r
}

# ============ 8) Heatmap integrat pe log2FC ============
# construim o matrice log2FC pentru toate contrastele (descriptiv)
log2fc_mat <- list()
for (nm in names(CONTRASTS)) {
  grp1 <- CONTRASTS[[nm]][1]; grp2 <- CONTRASTS[[nm]][2]
  resf <- file.path(DIR_DE, paste0(nm, "_", grp1, "_vs_", grp2, ".csv"))
  if (file.exists(resf)) {
    df <- readr::read_csv(resf, show_col_types = FALSE)
    log2fc_mat[[nm]] <- df %>% select(miRNA, log2FC) %>% tibble::column_to_rownames("miRNA")
  }
}
if (length(log2fc_mat) > 0) {
  L <- Reduce(function(x,y) { xn <- union(rownames(x), rownames(y));
                              x2 <- x[match(xn, rownames(x)), , drop=FALSE]; rownames(x2) <- xn
                              y2 <- y[match(xn, rownames(y)), , drop=FALSE]; rownames(y2) <- xn
                              cbind(x2, y2) }, log2fc_mat)
  colnames(L) <- names(log2fc_mat)
  L[is.na(L)] <- 0
  write_csv(as.data.frame(L), file.path(DIR_DE, "log2FC_all_contrasts.csv"))

  # top după varianță
  var_idx <- order(matrixStats::rowVars(as.matrix(L)), decreasing = TRUE)[seq_len(min(TOP_N_VAR, nrow(L)))]
  L_top <- L[var_idx, , drop=FALSE]
  png(file.path(DIR_PLOTS, "global_heatmap_top30.png"), width=1400, height=1200, res=150)
  pheatmap::pheatmap(L_top, main = "Global log2FC (top var miRNA)", color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(255))
  dev.off()

  # praguri utile (1.5x, 2x, 2.5x)
  for (thr in LOG2FC_THRESHES) {
    cnts <- colSums(abs(L) >= thr)
    write_csv(data.frame(contrast = names(cnts), n = as.integer(cnts)),
              file.path(DIR_DE, sprintf("counts_log2FC_ge_%.3f.csv", thr)))
  }
}

message("=== DONE ===")
message("Outputs written under: ", OUT_DIR)
