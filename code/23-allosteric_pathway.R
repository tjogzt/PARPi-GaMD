#!/usr/bin/env Rscript
# 23-allosteric_pathway.R — Allosteric pathway analysis & difference contact maps
# 
# 1A: Shortest-path analysis on DCCM networks from binding pocket → HD-ART interface
#      Compare AZD5305 vs Type II inhibitors
# 1B: Difference contact maps AZD5305 vs talazoparib/veliparib
#
# Input:  results/analysis/sys2_*_dccm.rds, ref_ca.pdb, CA DCDs
# Output: results/figures/Fig_Allosteric_Pathway.pdf
#         results/figures/Fig_Contact_Diff.pdf

library(bio3d)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
dcd_dir   <- "/Volumes/tjogzt4T/PARPi_data/md_analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Domain boundaries (consistent with 20-dccm_network.R)
DOMAINS <- list(
  ZnF1 = 1:96, ZnF2 = 97:214, ZnF3 = 215:383,
  BRCT = 384:517, WGR = 518:661,
  HD = 662:787, ART = 788:1014
)

INHIBITORS <- c(
  "sys2_APO" = "APO",
  "sys2_AZD5305" = "AZD5305",
  "sys2_talazoparib" = "Talazoparib",
  "sys2_veliparib" = "Veliparib",
  "sys2_niraparib" = "Niraparib",
  "sys2_olaparib" = "Olaparib",
  "sys2_rucaparib" = "Rucaparib"
)

TYPE <- c(
  "sys2_APO" = "APO",
  "sys2_AZD5305" = "Type_II",
  "sys2_talazoparib" = "Type_II",
  "sys2_niraparib" = "Type_II",
  "sys2_olaparib" = "Type_II",
  "sys2_rucaparib" = "Type_II",
  "sys2_veliparib" = "Type_III"
)

# ---- Theme -----------------------------------------------------------------
theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"))

# ---- Load reference PDB ----------------------------------------------------
ref_pdb <- read.pdb(file.path(data_dir, "ref_ca.pdb"))
ca_inds <- atom.select(ref_pdb, elety = "CA")
n_ca <- length(ca_inds$atom)
resno_vec <- ref_pdb$atom$resno[ca_inds$atom]

resid_to_domain <- function(r) {
  for (nm in names(DOMAINS)) {
    if (r %in% DOMAINS[[nm]]) return(nm)
  }
  return(NA_character_)
}
resid_dom <- sapply(resno_vec, resid_to_domain)

# Domain colors
dom_colors <- c(
  ZnF1 = "#66C2A5", ZnF2 = "#FC8D62", ZnF3 = "#8DA0CB",
  BRCT = "#E78AC3", WGR = "#A6D854",
  HD = "#FFD92F", ART = "#E5C494"
)

# ==============================================================================
# 1A: Shortest-path allosteric pathway analysis
# ==============================================================================

cat("========== 1A: Allosteric Pathway Analysis ==========\n\n")

# ---- Define CAT binding pocket residues (5 A from talazoparib in 4DQY) ----
pocket_resnos <- c(
  763, 764, 766, 767, 769, 770, 785, 786, 787, 788, 789, 790, 791,
  792, 793, 794, 795, 840, 841, 842, 843, 844, 845, 846, 847,
  869, 870, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881,
  907, 908
)

# Map to CA indices
pocket_ca_idx <- which(resno_vec %in% pocket_resnos)
cat(sprintf("Binding pocket: %d residues -> %d CA atoms\n",
            length(pocket_resnos), length(pocket_ca_idx)))

# ---- Define HD-ART interface residues (CA-CA < 8 A between domains) ----
hd_idx <- which(resid_dom == "HD")
art_idx <- which(resid_dom == "ART")
ref_xyz <- ref_pdb$xyz[ca_inds$xyz]

# Compute CA-CA distance matrix for HD vs ART
hd_art_dist <- matrix(NA, length(hd_idx), length(art_idx))
for (i in seq_along(hd_idx)) {
  for (j in seq_along(art_idx)) {
    xi <- ref_xyz[ (hd_idx[i]-1)*3 + 1:3 ]
    xj <- ref_xyz[ (art_idx[j]-1)*3 + 1:3 ]
    hd_art_dist[i, j] <- sqrt(sum((xi - xj)^2))
  }
}

interface_pairs <- which(hd_art_dist < 8.0, arr.ind = TRUE)
interface_hd <- hd_idx[unique(interface_pairs[, 1])]
interface_art <- art_idx[unique(interface_pairs[, 2])]
interface_idx <- sort(unique(c(interface_hd, interface_art)))
cat(sprintf("HD-ART interface: %d residues (CA-CA < 8 A)\n", length(interface_idx)))
cat(sprintf("  HD side: %d, ART side: %d\n",
            length(unique(interface_pairs[, 1])),
            length(unique(interface_pairs[, 2]))))

# ---- Build DCCM network and compute shortest paths ----
analyze_pathways <- function(sys_name) {
  dccm_file <- file.path(data_dir, paste0(sys_name, "_dccm.rds"))
  if (!file.exists(dccm_file)) {
    cat(sprintf("[SKIP] %s: no DCCM cache\n", sys_name))
    return(NULL)
  }
  
  cij <- readRDS(dccm_file)
  
  # Build network: edge weight = 1 - |correlation| (shorter path = stronger correlation)
  cij_abs <- abs(cij)
  diag(cij_abs) <- 0
  
  # Distance matrix (1 - |r|)
  dist_mat <- 1 - cij_abs
  dist_mat[dist_mat < 0] <- 0  # safety: avoid negative distances
  
  g <- graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Shortest paths from each pocket residue to each interface residue
  path_results <- list()
  path_count <- 0
  
  for (src in pocket_ca_idx) {
    sp <- shortest_paths(g, from = src, to = interface_idx, output = "both")
    
    for (k in seq_along(interface_idx)) {
      if (length(sp$vpath[[k]]) == 0) next  # no path
      
      path_count <- path_count + 1
      path_nodes <- as.numeric(sp$vpath[[k]])
      ep <- sp$epath[[k]]
      # Handle direct paths (path_length == 1: self-loop)
      if (length(ep) == 0) {
        path_dist <- 0
      } else {
        path_dist <- ep$weight
      }
      
      path_results[[path_count]] <- data.frame(
        source = resno_vec[src],
        target = resno_vec[interface_idx[k]],
        path_length = length(path_nodes),
        total_distance = sum(path_dist),
        mean_correlation = if (length(path_dist) > 0) 1 - mean(path_dist) else 1,
        n_hd = sum(resid_dom[path_nodes] == "HD", na.rm = TRUE),
        n_art = sum(resid_dom[path_nodes] == "ART", na.rm = TRUE),
        n_interface = sum(path_nodes %in% interface_idx),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(path_results) == 0) {
    cat(sprintf("  [WARN] %s: no valid paths found\n", sys_name))
    return(NULL)
  }
  
  paths_df <- bind_rows(path_results)
  
  # Summary statistics
  stats <- data.frame(
    system = sys_name,
    label = INHIBITORS[sys_name],
    type = TYPE[sys_name],
    n_paths = nrow(paths_df),
    mean_path_len = mean(paths_df$path_length),
    median_path_len = median(paths_df$path_length),
    mean_path_corr = mean(paths_df$mean_correlation),
    min_path_len = min(paths_df$path_length),
    n_direct = sum(paths_df$path_length == 2),  # direct edge
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("  %s: %d paths, mean length=%.1f, mean corr=%.3f, %d direct\n",
              INHIBITORS[sys_name], nrow(paths_df), 
              mean(paths_df$path_length), mean(paths_df$mean_correlation),
              sum(paths_df$path_length == 2)))
  
  list(paths = paths_df, stats = stats, graph = g)
}

pathway_results <- lapply(names(INHIBITORS), analyze_pathways)
names(pathway_results) <- names(INHIBITORS)

# Aggregate statistics
stats_all <- do.call(rbind, lapply(pathway_results, function(x) if (!is.null(x)) x$stats))
cat("\n--- Pathway Summary ---\n")
print(stats_all[, c("label", "type", "n_paths", "mean_path_len", "mean_path_corr")])

# ---- Identify unique AZD5305 pathways ----
azd <- pathway_results[["sys2_AZD5305"]]
tala <- pathway_results[["sys2_talazoparib"]]

if (!is.null(azd) && !is.null(tala)) {
  azd_paths <- azd$paths
  tala_paths <- tala$paths
  
  azd_by_src <- azd_paths %>% group_by(source) %>% 
    summarise(mean_len = mean(path_length), mean_corr = mean(mean_correlation), .groups = "drop")
  tala_by_src <- tala_paths %>% group_by(source) %>%
    summarise(mean_len = mean(path_length), mean_corr = mean(mean_correlation), .groups = "drop")
  
  merged <- merge(azd_by_src, tala_by_src, by = "source", suffixes = c("_azd", "_tala"))
  merged$delta_len <- merged$mean_len_azd - merged$mean_len_tala
  merged$delta_corr <- merged$mean_corr_azd - merged$mean_corr_tala
  
  cat(sprintf("\nAZD5305 vs Talazoparib pathway comparison:\n"))
  cat(sprintf("  delta_mean_path_len: %.2f (AZD: %.2f vs Tala: %.2f)\n",
              mean(merged$delta_len), mean(merged$mean_len_azd), mean(merged$mean_len_tala)))
  cat(sprintf("  delta_mean_corr: %.3f (AZD: %.3f vs Tala: %.3f)\n",
              mean(merged$delta_corr), mean(merged$mean_corr_azd), mean(merged$mean_corr_tala)))
  
  shorter_azd <- merged[merged$delta_len < -0.5, ]
  if (nrow(shorter_azd) > 0) {
    cat(sprintf("  %d pocket residues with shorter AZD5305 paths:\n", nrow(shorter_azd)))
    print(shorter_azd[, c("source", "mean_len_azd", "mean_len_tala", "delta_len")])
  }
}

# ---- Betweenness centrality hubs ----
btwn_all <- list()
for (sys_name in names(INHIBITORS)) {
  dccm_file <- file.path(data_dir, paste0(sys_name, "_dccm.rds"))
  if (!file.exists(dccm_file)) next
  cij <- readRDS(dccm_file)
  cij_abs <- abs(cij)
  diag(cij_abs) <- 0
  dist_mat <- 1 - cij_abs
  g <- graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  btwn_all[[sys_name]] <- betweenness(g, normalized = TRUE)
}

btwn_mat <- do.call(rbind, lapply(names(btwn_all), function(nm) btwn_all[[nm]]))
btwn_mean <- colMeans(btwn_mat)
top_hubs <- order(btwn_mean, decreasing = TRUE)[1:20]
cat("\n--- Top 20 Allosteric Hubs (mean betweenness across all systems) ---\n")
for (i in top_hubs) {
  cat(sprintf("  %d (%s): %.4f\n", resno_vec[i], resid_dom[i], btwn_mean[i]))
}

write.csv(stats_all, file.path(out_dir, "allosteric_pathway_stats.csv"), row.names = FALSE)
cat("\nSaved: allosteric_pathway_stats.csv\n")

# ---- Figure: Path length comparison ----
p_path_len <- ggplot(stats_all, aes(x = reorder(label, mean_path_len), 
                                      y = mean_path_len, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", mean_path_len)), 
            hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("APO" = "#999999", "Type_II" = "#2166AC", "Type_III" = "#B2182B")) +
  coord_flip(ylim = c(0, max(stats_all$mean_path_len) * 1.15)) +
  labs(x = NULL, y = "Mean Shortest Path Length",
       title = "Allosteric Path Length: Pocket -> HD-ART Interface",
       fill = NULL) +
  theme_7pt + theme(legend.position = "bottom")

p_path_corr <- ggplot(stats_all, aes(x = reorder(label, mean_path_corr),
                                       y = mean_path_corr, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_path_corr)), 
            hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("APO" = "#999999", "Type_II" = "#2166AC", "Type_III" = "#B2182B")) +
  coord_flip(ylim = c(0, max(stats_all$mean_path_corr) * 1.15)) +
  labs(x = NULL, y = "Mean Path Correlation",
       title = "Allosteric Coupling Strength: Pocket -> Interface",
       fill = NULL) +
  theme_7pt + theme(legend.position = "none")

# ==============================================================================
# 1B: Difference Contact Maps (from CA DCDs)
# ==============================================================================

cat("\n========== 1B: Difference Contact Maps ==========\n\n")

CONTACT_CUTOFF <- 8.0

compute_contact_map <- function(sys_name, max_frames = 200) {
  ca_dcd <- file.path(dcd_dir, paste0(sys_name, "_ca.dcd"))
  if (!file.exists(ca_dcd)) {
    cat(sprintf("[SKIP] %s: no CA DCD\n", sys_name))
    return(NULL)
  }
  
  cat(sprintf("  Loading %s...\n", ca_dcd))
  dcd <- read.dcd(ca_dcd)
  n_total <- dim(dcd)[1]
  n_atoms <- dim(dcd)[2] / 3
  
  # Subsample evenly
  if (n_total > max_frames) {
    stride <- max(floor(n_total / max_frames), 1)
    frame_idx <- seq(1, n_total, by = stride)
    if (length(frame_idx) > max_frames) frame_idx <- frame_idx[1:max_frames]
  } else {
    frame_idx <- seq_len(n_total)
  }
  n_frames <- length(frame_idx)
  cat(sprintf("  Subsampled %d/%d frames\n", n_frames, n_total))
  
  cmap <- matrix(0, n_atoms, n_atoms)
  
  for (t in frame_idx) {
    xyz_mat <- matrix(dcd[t, ], ncol = 3, byrow = TRUE)
    d_mat <- as.matrix(dist(xyz_mat))
    cmap <- cmap + (d_mat < CONTACT_CUTOFF) * 1
  }
  
  cmap <- cmap / n_frames
  diag(cmap) <- 0
  cat(sprintf("  Done: %d frames, mean contact freq = %.3f\n", n_frames, mean(cmap[upper.tri(cmap)])))
  cmap
}

hd_art_range <- which(resno_vec >= 662 & resno_vec <= 1011)
cat(sprintf("HD-ART CA region: %d atoms\n", length(hd_art_range)))

systems_to_map <- c("sys2_AZD5305", "sys2_talazoparib", "sys2_veliparib", "sys2_APO")
cmaps <- list()

for (sys in systems_to_map) {
  cmaps[[sys]] <- compute_contact_map(sys)
}

# Compute delta maps
delta_AZD_tala <- NULL
delta_AZD_veli <- NULL
delta_tala_veli <- NULL

if (!is.null(cmaps[["sys2_AZD5305"]]) && !is.null(cmaps[["sys2_talazoparib"]])) {
  delta_AZD_tala <- cmaps[["sys2_AZD5305"]] - cmaps[["sys2_talazoparib"]]
}
if (!is.null(cmaps[["sys2_AZD5305"]]) && !is.null(cmaps[["sys2_veliparib"]])) {
  delta_AZD_veli <- cmaps[["sys2_AZD5305"]] - cmaps[["sys2_veliparib"]]
}
if (!is.null(cmaps[["sys2_talazoparib"]]) && !is.null(cmaps[["sys2_veliparib"]])) {
  delta_tala_veli <- cmaps[["sys2_talazoparib"]] - cmaps[["sys2_veliparib"]]
}

plot_delta_contact <- function(delta_mat, title_str) {
  delta_sub <- delta_mat[hd_art_range, hd_art_range]
  resnos_sub <- resno_vec[hd_art_range]
  
  delta_df <- expand.grid(i = seq_len(nrow(delta_sub)), j = seq_len(ncol(delta_sub)))
  delta_df$res_i <- resnos_sub[delta_df$i]
  delta_df$res_j <- resnos_sub[delta_df$j]
  delta_df$delta <- as.vector(delta_sub)
  delta_df <- delta_df[delta_df$i > delta_df$j, ]
  
  max_abs <- max(abs(delta_df$delta), 0.05)
  
  ggplot(delta_df, aes(x = res_i, y = res_j, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-max_abs, max_abs),
                         name = expression(Delta*"Contact\nFreq"),
                         oob = scales::squish) +
    geom_vline(xintercept = 787.5, color = "grey40", linewidth = 0.3, linetype = "dashed") +
    geom_hline(yintercept = 787.5, color = "grey40", linewidth = 0.3, linetype = "dashed") +
    annotate("text", x = 725, y = 995, label = "HD", size = 2.5, color = "grey30") +
    annotate("text", x = 900, y = 995, label = "ART", size = 2.5, color = "grey30") +
    annotate("text", x = 725, y = 700, label = "HD", size = 2.5, color = "grey30") +
    annotate("text", x = 900, y = 700, label = "ART", size = 2.5, color = "grey30") +
    coord_fixed() +
    labs(x = "Residue", y = "Residue", title = title_str) +
    theme_7pt + theme(axis.text = element_text(size = 4),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

p_delta_list <- list()

if (!is.null(delta_AZD_tala)) {
  p_delta_list[["AZD_Tala"]] <- plot_delta_contact(delta_AZD_tala,
    "Delta Contact: AZD5305 - Talazoparib")
}
if (!is.null(delta_AZD_veli)) {
  p_delta_list[["AZD_Veli"]] <- plot_delta_contact(delta_AZD_veli,
    "Delta Contact: AZD5305 - Veliparib")
}
if (!is.null(delta_tala_veli)) {
  p_delta_list[["Tala_Veli"]] <- plot_delta_contact(delta_tala_veli,
    "Delta Contact: Talazoparib - Veliparib")
}

# Export significant differential contacts
if (!is.null(delta_AZD_tala)) {
  delta_hd_art <- delta_AZD_tala[hd_art_range, hd_art_range]
  resnos_sub <- resno_vec[hd_art_range]
  up_idx <- upper.tri(delta_hd_art)
  
  sig_df <- data.frame(
    res_i = resnos_sub[row(delta_hd_art)[up_idx]],
    res_j = resnos_sub[col(delta_hd_art)[up_idx]],
    dom_i = resid_dom[hd_art_range][row(delta_hd_art)[up_idx]],
    dom_j = resid_dom[hd_art_range][col(delta_hd_art)[up_idx]],
    delta_AZD_Tala = delta_hd_art[up_idx],
    AZD5305_freq = cmaps[["sys2_AZD5305"]][hd_art_range, hd_art_range][up_idx],
    Talazoparib_freq = cmaps[["sys2_talazoparib"]][hd_art_range, hd_art_range][up_idx]
  )
  
  sig_df <- sig_df[order(-abs(sig_df$delta_AZD_Tala)), ]
  write.csv(sig_df, file.path(out_dir, "AZD5305_Talazoparib_contact_diff.csv"), row.names = FALSE)
  cat("Saved: AZD5305_Talazoparib_contact_diff.csv\n")
  
  top_hits <- head(sig_df, 10)
  cat("\nTop 10 differential contacts AZD5305 vs Talazoparib:\n")
  for (k in seq_len(nrow(top_hits))) {
    cat(sprintf("  %d(%s) <-> %d(%s): delta=%.3f (AZD=%.3f, Tala=%.3f)\n",
                top_hits$res_i[k], top_hits$dom_i[k],
                top_hits$res_j[k], top_hits$dom_j[k],
                top_hits$delta_AZD_Tala[k],
                top_hits$AZD5305_freq[k],
                top_hits$Talazoparib_freq[k]))
  }
}

# ---- Figures ----
fig_pathway <- (p_path_len / p_path_corr) +
  plot_annotation(
    title = "Allosteric Pathway Analysis: Binding Pocket to HD-ART Interface",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

cairo_pdf(file.path(out_dir, "Fig_Allosteric_Pathway.pdf"),
          width = 190/25.4, height = 180/25.4, pointsize = 7)
print(fig_pathway)
dev.off()
cat("Saved: Fig_Allosteric_Pathway.pdf\n")

if (length(p_delta_list) > 0) {
  fig_contact <- wrap_plots(p_delta_list, ncol = 1) +
    plot_annotation(
      title = "Residue Contact Frequency Difference Maps (CA-CA < 8 A)",
      theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
  
  cairo_pdf(file.path(out_dir, "Fig_Contact_Diff.pdf"),
            width = 190/25.4, height = 200/25.4, pointsize = 7)
  print(fig_contact)
  dev.off()
  cat("Saved: Fig_Contact_Diff.pdf\n")
}

message("===== Allosteric pathway + contact map analysis complete =====")
