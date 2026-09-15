#!/usr/bin/env Rscript
# 20-dccm_network.R — DCCM + community network analysis for S2 systems
# Computes residue-residue cross-correlation and identifies
# allosteric signal transduction pathways (Fig 4D)
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

# Domain boundaries
DOMAINS <- list(
  ZnF1 = 1:96, ZnF2 = 97:214, ZnF3 = 215:383,
  BRCT = 384:517, WGR = 518:661,
  HD = 662:787, ART = 788:1014
)

# Inhibitor classification
INHIBITORS <- data.frame(
  system = c("sys2_APO", "sys2_AZD5305", "sys2_talazoparib",
             "sys2_veliparib", "sys2_niraparib", "sys2_olaparib", "sys2_rucaparib"),
  label  = c("APO", "AZD5305", "Talazoparib",
             "Veliparib", "Niraparib", "Olaparib", "Rucaparib"),
  type   = c("APO", "Type II", "Type II",
             "Type III", "Type III", "Type III", "Type III"),
  stringsAsFactors = FALSE
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
cat(sprintf("Reference: %d CA atoms\n", n_ca))

# ---- Domain color mapping --------------------------------------------------
dom_colors <- c(
  ZnF1 = "#66C2A5", ZnF2 = "#FC8D62", ZnF3 = "#8DA0CB",
  BRCT = "#E78AC3", WGR = "#A6D854",
  HD = "#FFD92F", ART = "#E5C494"
)

resid_to_domain <- function(r) {
  for (nm in names(DOMAINS)) {
    if (r %in% DOMAINS[[nm]]) return(nm)
  }
  return(NA_character_)
}
resno_vec <- ref_pdb$atom$resno[ca_inds$atom]
resid_dom <- unlist(sapply(resno_vec, resid_to_domain))
resid_col <- dom_colors[resid_dom]

# ---- Process each system ---------------------------------------------------
compute_dccm <- function(sys_name, ca_dcd_path) {
  cat(sprintf("\n--- %s ---\n", sys_name))
  
  # Load or compute DCCM
  dccm_file <- file.path(data_dir, paste0(sys_name, "_dccm.rds"))
  if (file.exists(dccm_file)) {
    cij <- readRDS(dccm_file)
    cat(sprintf("  Loaded cached: %d x %d\n", nrow(cij), ncol(cij)))
  } else {
    dcd <- read.dcd(ca_dcd_path)
    n_frames <- dim(dcd)[1]
    cat(sprintf("  Frames: %d\n", n_frames))
    cij <- dccm(dcd)
    saveRDS(cij, dccm_file)
    cat(sprintf("  Saved: %s\n", dccm_file))
  }
  
  # Community network analysis
  # Convert DCCM to network (edges = |correlation| > cutoff)
  cutoff <- 0.4
  cij_abs <- abs(cij)
  diag(cij_abs) <- 0
  
  # Create weighted adjacency matrix
  adj <- cij_abs
  adj[adj < cutoff] <- 0
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Community detection (Louvain)
  comm <- cluster_louvain(g)
  n_comm <- length(unique(comm$membership))
  cat(sprintf("  Communities: %d (cutoff=%.1f)\n", n_comm, cutoff))
  
  # Modularity
  mod <- modularity(comm)
  cat(sprintf("  Modularity: %.3f\n", mod))
  
  # Betweenness centrality for key residues
  btwn <- betweenness(g, normalized = TRUE)
  top_btwn <- order(btwn, decreasing = TRUE)[1:10]
  cat("  Top 10 betweenness residues:\n")
  for (i in top_btwn) {
    dom <- resid_dom[i]
    cat(sprintf("    %d (%s): %.4f\n", resno_vec[i], dom, btwn[i]))
  }
  
  # Cross-domain correlation summary
  hd_idx <- which(resid_dom == "HD")
  art_idx <- which(resid_dom == "ART")
  hd_art_corr <- mean(cij_abs[hd_idx, art_idx])
  cat(sprintf("  Mean |corr| HD↔ART: %.4f\n", hd_art_corr))
  
  # Return results
  list(
    system = sys_name,
    cij = cij,
    graph = g,
    communities = comm,
    modularity = mod,
    betweenness = btwn,
    hd_art_corr = hd_art_corr,
    n_comm = n_comm
  )
}

# Process available systems
results <- list()
summary_rows <- list()

for (i in seq_len(nrow(INHIBITORS))) {
  sys <- INHIBITORS$system[i]
  ca_dcd <- file.path(dcd_dir, paste0(sys, "_ca.dcd"))
  
  if (!file.exists(ca_dcd)) {
    cat(sprintf("[SKIP] %s: no CA DCD\n", sys))
    next
  }
  
  res <- compute_dccm(sys, ca_dcd)
  results[[sys]] <- res
  
  summary_rows[[sys]] <- data.frame(
    system = sys,
    label = INHIBITORS$label[i],
    type = INHIBITORS$type[i],
    modularity = res$modularity,
    n_communities = res$n_comm,
    hd_art_corr = res$hd_art_corr,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
cat("\n========== DCCM SUMMARY ==========\n")
print(summary_df)

# ---- Generate DCCM heatmap figure (Type II vs Type III) --------------------
# Compare talazoparib (Type II) vs veliparib (Type III)
tala <- results[["sys2_talazoparib"]]
veli <- results[["sys2_veliparib"]]

if (!is.null(tala) && !is.null(veli)) {
  # Delta DCCM
  delta_dccm <- tala$cij - veli$cij
  
  # Focus on HD-ART region
  hd_art_range <- which(resid_dom %in% c("HD", "ART"))
  delta_sub <- delta_dccm[hd_art_range, hd_art_range]
  resnos_sub <- resno_vec[hd_art_range]
  
  # Convert to long format for ggplot
  delta_df <- expand.grid(i = seq_len(nrow(delta_sub)), j = seq_len(ncol(delta_sub)))
  delta_df$res_i <- resnos_sub[delta_df$i]
  delta_df$res_j <- resnos_sub[delta_df$j]
  delta_df$delta <- as.vector(delta_sub)
  delta_df <- delta_df[delta_df$i < delta_df$j, ]  # upper triangle only
  
  # Domain annotation for axes
  delta_df$dom_i <- resid_dom[hd_art_range][delta_df$i]
  delta_df$dom_j <- resid_dom[hd_art_range][delta_df$j]
  
  p_delta <- ggplot(delta_df, aes(x = res_i, y = res_j, fill = delta)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = expression(Delta*"DCCM"),
                         limits = c(-0.5, 0.5), oob = scales::squish) +
    # Domain boundary lines
    geom_vline(xintercept = 787.5, color = "grey40", linewidth = 0.3, linetype = "dashed") +
    geom_hline(yintercept = 787.5, color = "grey40", linewidth = 0.3, linetype = "dashed") +
    annotate("text", x = 725, y = 975, label = "HD", size = 2.5, color = "grey30") +
    annotate("text", x = 900, y = 975, label = "ART", size = 2.5, color = "grey30") +
    annotate("text", x = 725, y = 730, label = "HD", size = 2.5, color = "grey30") +
    annotate("text", x = 900, y = 730, label = "ART", size = 2.5, color = "grey30") +
    coord_fixed() +
    labs(x = "Residue", y = "Residue",
         title = "ΔDCCM: Talazoparib - Veliparib (HD-ART)") +
    theme_7pt + theme(axis.text = element_text(size = 4),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # ---- Modularity bar plot ----
  p_mod <- ggplot(summary_df, aes(x = reorder(label, modularity), y = modularity, fill = type)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("APO" = "#666666", "Type II" = "#2166AC", "Type III" = "#B2182B")) +
    labs(x = NULL, y = "Modularity", title = "Community Modularity by Inhibitor",
         fill = NULL) +
    theme_7pt + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = "bottom")
  
  # ---- HD-ART correlation bar plot ----
  p_corr <- ggplot(summary_df, aes(x = reorder(label, hd_art_corr), y = hd_art_corr, fill = type)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("APO" = "#666666", "Type II" = "#2166AC", "Type III" = "#B2182B")) +
    labs(x = NULL, y = "Mean |corr|", title = "HD↔ART Cross-Correlation",
         fill = NULL) +
    theme_7pt + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = "none")
  
  # ---- Assemble ----
  master <- (p_delta | (p_mod / p_corr)) +
    plot_layout(widths = c(1.5, 0.8)) +
    plot_annotation(
      title = "Dynamic Network Analysis: Talazoparib Strengthens HD-ART Coupling",
      theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
  
  cairo_pdf(file.path(out_dir, "Fig_DCCM_Network.pdf"),
            width = 190/25.4, height = 150/25.4, pointsize = 7)
  print(master)
  dev.off()
  cat("\nSaved: Fig_DCCM_Network.pdf\n")
}

# ---- Key findings ----
cat("\n========== KEY DCCM FINDINGS ==========\n")
cat(sprintf("Type II mean modularity: %.3f\n",
            mean(summary_df$modularity[summary_df$type == "Type II"])))
cat(sprintf("Type III mean modularity: %.3f\n",
            mean(summary_df$modularity[summary_df$type == "Type III"])))
cat(sprintf("Type II HD-ART correlation: %.3f\n",
            mean(summary_df$hd_art_corr[summary_df$type == "Type II"])))
cat(sprintf("Type III HD-ART correlation: %.3f\n",
            mean(summary_df$hd_art_corr[summary_df$type == "Type III"])))

message("===== DCCM analysis complete =====")
