# ===============================================================
# META-ANALYSIS RESULT VISUALIZATION (GGSEG)
# ---------------------------------------------------------------
# This script visualizes ROI-level meta-analysis results on cortical
# and subcortical atlases using ggseg. For each session (RS, WM, FI),
# it:
#   - computes FDR-corrected p-values by session,
#   - maps pFDR and meta z-values (zval) to ggseg atlases,
#   - exports multi-panel PDFs.
#
# Optional:
#   - Δ (delta) & heterogeneity (QE) maps,
#   - cross-diagnosis scatter comparisons (SCZ vs BD).
#
# Author: Roberta Passiatore
# Updated: 2025-10-06
# ===============================================================

# ---- Packages ----
need <- c("ggseg", "ggsegSchaefer", "ggplot2", "dplyr", "patchwork", "grid")
to_install <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(ggseg)
library(ggsegSchaefer)
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

# ---- CONFIG ----
# Cortical inputs
cort_results_file <- "res.GroupDiff_NCvsBD_corticalROI_REML_standard_metanalysis.csv"
cort_roi_names    <- "/Users/roberta.passiatore/Documents/Motion-Profiling/analyses/roi-names.csv"
cort_atlas_name   <- "schaefer17_100"
cort_output_pdf   <- "Metanalysis-results_corticalROI_NCvsBD_ME_pFDR_z.pdf"

# Subcortical inputs
sub_results_file  <- "res.GroupDiff_NCvsBD_subcorticalROI_REML_standard_metanalysis.csv"
sub_roi_names     <- "/Users/roberta.passiatore/Documents/Motion-Profiling/analyses/roi-names-sub.csv"
sub_atlas_name    <- "aseg"
sub_output_pdf    <- "Metanalysis-results_subcorticalROI_NCvsBD_ME_pFDR_z.pdf"

# Optional Delta/QE maps
do_delta_qe_cort  <- TRUE
delta_cort_file   <- "res.GroupDiff_NCvsBD_corticalROI_REML_delta_metanalysis.csv"
delta_cort_pdf    <- "Metanalysis-results_corticalROI_NCvsBD_ME_Delta_QE.pdf"

do_delta_qe_sub   <- TRUE
delta_sub_file    <- "res.GroupDiff_NCvsBD_subcorticalROI_REML_delta_metanalysis.csv"
delta_sub_pdf     <- "Metanalysis-results_subcorticalROI_NCvsBD_ME_Delta_QE.pdf"

# Cross-diagnosis scatter (SCZ~BD)
do_diag_compare   <- TRUE
bd_cort_delta_file <- "res.NCvsBD_GroupDiff_NC_corticalROI_delta_metanalysis.csv"
bd_sub_delta_file  <- "res.NCvsBD_GroupDiff_NC_subcorticalROI_delta_metanalysis.csv"
scz_cort_delta_file <- "res.GroupDiff_NCvsSCZ_corticalROI_delta_metanalysis.csv"
scz_sub_delta_file  <- "res.GroupDiff_NCvsSCZ_subcorticalROI_delta_metanalysis.csv"

# FDR alpha and plotting ranges
alpha_fdr <- 0.05
z_limits_cort <- c(-15, 15)
z_limits_sub  <- c(-10, 10)
qe_limits_cort <- c(0, 20)
qe_limits_sub  <- c(0, 10)

# Sessions to display (must match your results)
sessions_order <- c("RS","WM","FI")

# ---- Helpers ----
session_fdr <- function(df, pcol = "pval", sess_col = "session") {
  # FDR per session
  df <- df %>% mutate(pFDR = NA_real_)
  for (s in unique(df[[sess_col]])) {
    idx <- which(df[[sess_col]] == s & is.finite(df[[pcol]]))
    if (length(idx)) df$pFDR[idx] <- p.adjust(df[[pcol]][idx], method = "fdr")
  }
  df
}

atlas_merge <- function(atlas_obj, map_df, key = "label") {
  at <- atlas_obj
  at$data <- dplyr::left_join(at$data, map_df, by = key)
  at
}

make_y_annotation <- function(text) {
  ylab <- textGrob(text, rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
  wrap_elements(full = ylab)
}

panel_grid3 <- function(p_RS, p_WM, p_FI, y_annot) {
  y_annot + (p_RS / p_WM / p_FI) +
    plot_layout(guides = "collect", widths = c(1, 10)) +
    plot_annotation(title = "")
}

# ---- Load atlases ----
data(list = cort_atlas_name, package = "ggsegSchaefer")
cort_atlas <- get(cort_atlas_name)
sub_atlas  <- get(sub_atlas_name)

# =========================================================
# CORTICAL: pFDR & z maps
# =========================================================
rois <- read.csv(cort_roi_names, stringsAsFactors = FALSE)
res  <- read.csv(cort_results_file, stringsAsFactors = FALSE)

# Attach atlas label by ROI
res$label <- rois$Atlas[match(res$ROI, rois$ROI)]

# FDR per session
res <- session_fdr(res, pcol = "pval", sess_col = "session")

# Build plotting DF and merge
plot_df <- res[, c("label","pFDR","zval","session")]
# mask pFDR > alpha
plot_df$zval[plot_df$pFDR > alpha_fdr] <- NA
plot_df$pFDR[plot_df$pFDR > alpha_fdr] <- NA

# One atlas per session
mk_cort_plots <- function(sess) {
  tmp <- plot_df %>% filter(session == sess)
  a1 <- atlas_merge(cort_atlas, tmp[, c("label","pFDR")])
  a2 <- atlas_merge(cort_atlas, tmp[, c("label","zval")])

  p_pFDR <- ggseg(atlas = a1, mapping = aes(fill = pFDR),
                  colour = "grey37", size = 0.3, position = "dispersed") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = "grey", limits = c(0, alpha_fdr)) +
    theme_brain() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 10)) +
    ggtitle(paste("pFDR", sess))

  p_z <- ggseg(atlas = a2, mapping = aes(fill = zval),
               colour = "grey37", size = 0.3, position = "dispersed") +
    scale_fill_gradient(low = "blue", high = "red", na.value = "grey", limits = z_limits_cort) +
    theme_brain() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 10)) +
    ggtitle(paste("z", sess))

  list(p_pFDR = p_pFDR, p_z = p_z)
}

plots_RS <- mk_cort_plots("RS")
plots_WM <- mk_cort_plots("WM")
plots_FI <- mk_cort_plots("FI")

ylab <- make_y_annotation("         Corrected        Standard")

combined_pFDR <- panel_grid3(plots_RS$p_pFDR, plots_WM$p_pFDR, plots_FI$p_pFDR, ylab)
combined_z    <- panel_grid3(plots_RS$p_z,    plots_WM$p_z,    plots_FI$p_z,    ylab)

pdf(file = cort_output_pdf, width = 10, height = 7)
print(combined_pFDR)
print(combined_z)
dev.off()

# =========================================================
# SUBCORTICAL: pFDR & z maps
# =========================================================
rois_sub <- read.csv(sub_roi_names, stringsAsFactors = FALSE)
res_sub  <- read.csv(sub_results_file, stringsAsFactors = FALSE)

# Add labels and mirror left/right if desired (optional trick shown originally)
res_sub$label <- rois_sub$Atlas[match(res_sub$ROI, rois_sub$ROI)]
# Optionally duplicate left/right to fill atlas symmetry
tmp_mirror <- res_sub
tmp_mirror$label <- gsub("Right", "Left", tmp_mirror$label)
res_sub <- rbind(res_sub, tmp_mirror)

res_sub <- session_fdr(res_sub, pcol = "pval", sess_col = "session")

plot_df_sub <- res_sub[, c("label","pFDR","zval","session")]
plot_df_sub$zval[plot_df_sub$pFDR > alpha_fdr] <- NA
plot_df_sub$pFDR[plot_df_sub$pFDR > alpha_fdr] <- NA

mk_sub_plots <- function(sess) {
  tmp <- plot_df_sub %>% filter(session == sess)
  a1 <- atlas_merge(sub_atlas, tmp[, c("label","pFDR")])
  a2 <- atlas_merge(sub_atlas, tmp[, c("label","zval")])

  p_pFDR <- ggseg(atlas = a1, mapping = aes(fill = pFDR),
                  colour = "grey37", size = 0.3, position = "dispersed") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = "grey", limits = c(0, alpha_fdr)) +
    theme_brain() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 10)) +
    ggtitle(paste("pFDR", sess))

  p_z <- ggseg(atlas = a2, mapping = aes(fill = zval),
               colour = "grey37", size = 0.3, position = "dispersed") +
    scale_fill_gradient(low = "blue", high = "red", na.value = "grey", limits = z_limits_sub) +
    theme_brain() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 10)) +
    ggtitle(paste("z", sess))

  list(p_pFDR = p_pFDR, p_z = p_z)
}

s_RS <- mk_sub_plots("RS")
s_WM <- mk_sub_plots("WM")
s_FI <- mk_sub_plots("FI")

ylab_sub <- make_y_annotation("   Corrected        Standard")

combined_pFDR_sub <- panel_grid3(s_RS$p_pFDR, s_WM$p_pFDR, s_FI$p_pFDR, ylab_sub)
combined_z_sub    <- panel_grid3(s_RS$p_z,    s_WM$p_z,    s_FI$p_z,    ylab_sub)

pdf(file = sub_output_pdf, width = 8, height = 6)
print(combined_pFDR_sub)
print(combined_z_sub)
dev.off()

# =========================================================
# OPTIONAL: CORTICAL Δ & QE maps
# =========================================================
if (do_delta_qe_cort) {
  res_d <- read.csv(delta_cort_file, stringsAsFactors = FALSE)
  res_d$label <- rois$Atlas[match(res_d$ROI, rois$ROI)]
  res_d <- session_fdr(res_d, pcol = "pval", sess_col = "session")

  plot_df <- res_d[, c("label","QEp","pFDR","zval","session")]
  plot_df$QEpFDR <- ave(plot_df$QEp, plot_df$session, FUN = function(p) p.adjust(p, "fdr"))
  plot_df$zval[plot_df$pFDR > alpha_fdr | is.na(plot_df$pFDR)] <- NA
  plot_df$QE <- NA_real_
  # If QE is present in your file, replace the NA fill:
  if ("QE" %in% names(res_d)) {
    plot_df$QE <- res_d$QE[match(plot_df$label, res_d$label)]
  }

  mk_cort_delta <- function(sess) {
    tmp <- plot_df %>% filter(session == sess)
    a_qe <- atlas_merge(cort_atlas, tmp[, c("label","QE")])
    a_d  <- atlas_merge(cort_atlas, tmp[, c("label","zval")])

    p_qe <- ggseg(atlas = a_qe, mapping = aes(fill = QE),
                  colour = "grey37", size = 0.3, position = "dispersed") +
      scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey", limits = qe_limits_cort) +
      theme_minimal() + ggtitle(paste("QE", sess))

    p_delta <- ggseg(atlas = a_d, mapping = aes(fill = zval),
                     colour = "grey37", size = 0.3, position = "dispersed") +
      scale_fill_gradient(low = "white", high = "darkgreen", na.value = "grey", limits = c(0, 6)) +
      theme_brain() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 10)) +
      labs(fill = expression(Delta)) + ggtitle(paste("Δ", sess))

    list(p_qe = p_qe, p_delta = p_delta)
  }

  d_RS <- mk_cort_delta("RS")
  d_WM <- mk_cort_delta("WM")
  d_FI <- mk_cort_delta("FI")

  ylab_delta <- make_y_annotation("    Delta")
  comb_qe <- panel_grid3(d_RS$p_qe, d_WM$p_qe, d_FI$p_qe, ylab_delta)
  comb_d  <- panel_grid3(d_RS$p_delta, d_WM$p_delta, d_FI$p_delta, ylab_delta)

  pdf(file = delta_cort_pdf, width = 10, height = 7)
  print(comb_qe); print(comb_d)
  dev.off()
}

# =========================================================
# OPTIONAL: SUBCORTICAL Δ & QE maps
# =========================================================
if (do_delta_qe_sub) {
  res_d <- read.csv(delta_sub_file, stringsAsFactors = FALSE)
  res_d$label <- rois_sub$Atlas[match(res_d$ROI, rois_sub$ROI)]
  # mirror L/R as before
  tmp_m <- res_d; tmp_m$label <- gsub("Right", "Left", tmp_m$label)
  res_d <- rbind(res_d, tmp_m)

  res_d <- session_fdr(res_d, pcol = "pval", sess_col = "session")

  plot_df <- res_d[, c("label","QEp","pFDR","zval","session")]
  plot_df$QEpFDR <- ave(plot_df$QEp, plot_df$session, FUN = function(p) p.adjust(p, "fdr"))
  plot_df$zval[plot_df$pFDR > alpha_fdr | is.na(plot_df$pFDR)] <- NA
  plot_df$QE <- if ("QE" %in% names(res_d)) res_d$QE else NA_real_

  mk_sub_delta <- function(sess) {
    tmp <- plot_df %>% filter(session == sess)
    a_qe <- atlas_merge(sub_atlas, tmp[, c("label","QE")])
    a_d  <- atlas_merge(sub_atlas, tmp[, c("label","zval")])

    p_qe <- ggseg(atlas = a_qe, mapping = aes(fill = QE),
                  colour = "grey37", size = 0.3, position = "dispersed") +
      scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey", limits = qe_limits_sub) +
      theme_minimal() + ggtitle(paste("QE", sess))

    p_delta <- ggseg(atlas = a_d, mapping = aes(fill = zval),
                     colour = "grey37", size = 0.3, position = "dispersed") +
      scale_fill_gradient(low = "white", high = "darkgreen", na.value = "grey", limits = c(0, 6)) +
      theme_brain() + theme(text = element_text(size = 16), axis.text.x = element_text(size = 10)) +
      labs(fill = expression(Delta)) + ggtitle(paste("Δ", sess))

    list(p_qe = p_qe, p_delta = p_delta)
  }

  d_RS <- mk_sub_delta("RS")
  d_WM <- mk_sub_delta("WM")
  d_FI <- mk_sub_delta("FI")

  ylab_delta <- make_y_annotation("    Delta")
  comb_qe <- panel_grid3(d_RS$p_qe, d_WM$p_qe, d_FI$p_qe, ylab_delta)
  comb_d  <- panel_grid3(d_RS$p_delta, d_WM$p_delta, d_FI$p_delta, ylab_delta)

  pdf(file = delta_sub_pdf, width = 8, height = 6)
  print(comb_qe); print(comb_d)
  dev.off()
}

# =========================================================
# OPTIONAL: SCZ ~ BD cross-diagnosis comparisons
# =========================================================
if (do_diag_compare) {
  suppressPackageStartupMessages({
    if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
    if (!requireNamespace("ggthemes", quietly = TRUE)) install.packages("ggthemes")
    library(ggpubr); library(ggthemes)
  })

  # Load and combine
  bd <- read.csv(bd_cort_delta_file, stringsAsFactors = FALSE) %>%
    rbind(read.csv(bd_sub_delta_file, stringsAsFactors = FALSE))
  scz <- read.csv(scz_cort_delta_file, stringsAsFactors = FALSE) %>%
    rbind(read.csv(scz_sub_delta_file, stringsAsFactors = FALSE))

  # FDR per session
  for (S in unique(bd$session)) {
    idx <- which(bd$session == S)
    if (length(idx)) bd$pFDR[idx] <- p.adjust(bd$pval[idx], "fdr")
  }
  for (S in unique(scz$session)) {
    idx <- which(scz$session == S)
    if (length(idx)) scz$pFDR[idx] <- p.adjust(scz$pval[idx], "fdr")
  }

  bd$Diagnosis <- "BD"; scz$Diagnosis <- "SCZ"
  names(bd)  <- paste0(names(bd), "_BD")
  names(scz) <- paste0(names(scz), "_SCZ")

  df <- cbind(scz, bd)  # aligns by row order; ensure files are aligned by ROI ordering upstream

  p1 <- ggplot(df, aes(x = Delta_SCZ, y = Delta_BD)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgreen") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4, parse = FALSE) +
    theme_few(base_size = 14) + xlim(-10, 10) + ylim(-10, 10) +
    labs(title = "SCZ ~ BD (Δ across ROIs)", x = "Δ z-score (SCZ)", y = "Δ z-score (BD)")

  p2 <- ggplot(df, aes(x = zval_SCZ, y = zval_BD)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgreen") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4, parse = FALSE) +
    theme_few(base_size = 14) + xlim(-10, 10) + ylim(-10, 10) +
    labs(title = "SCZ ~ BD (meta z across ROIs)", x = "Meta z (SCZ)", y = "Meta z (BD)")

  print(p1 + p2)
}
