#!/usr/bin/env Rscript
# 22-cat_pocket_mapping.R — CAT pocket residue mapping (Fig 5A/5B)
# Maps differential RMSF between Type II and III onto CAT pocket residues
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# CAT pocket residues from PDB 4DQY — inhibitor binding site
# ART catalytic pocket: residues within 5Å of talazoparib in 4DQY
CAT_POCKET <- c(
  # NAD+ binding / catalytic triad
  859, 860, 861, 862, 863,  # G-rich loop
  872, 873, 874,             # catalytic loop
  886, 887, 888, 889, 890,   # D-loop (D887 catalytic)
  896, 897, 898,             # Y896 (stacking)
  903, 904, 905, 906, 907, 908,  # DON loop
  # Pocket floor/walls
  788, 789, 790, 791,        # N-terminal ART
  840, 841, 842, 843,        # helix
  850, 851, 852, 853,        # loop
  866, 867, 868, 869, 870, 871  # beta-sheet
)

CAT_FEATURES <- data.frame(
  resid = c(862, 863, 872, 887, 888, 890, 896, 903, 904, 907, 908, 788),
  label = c("G862", "R863", "S872", "D887", "S888", "D890", "Y896", "Y903", "S904", "Y907", "E908", "L788"),
  feature = c("G-loop", "G-loop", "Cat loop", "Cat D", "D-loop", "D-loop",
              "Stacking", "DON", "DON", "DON", "DON", "N-term"),
  stringsAsFactors = FALSE
)

# ---- Load RMSF data --------------------------------------------------------
TYPE_MAP <- data.frame(
  system = c("sys2_APO", "sys2_AZD5305", "sys2_talazoparib",
             "sys2_veliparib", "sys2_niraparib", "sys2_olaparib", "sys2_rucaparib"),
  label  = c("APO", "AZD5305", "Talazoparib",
             "Veliparib", "Niraparib", "Olaparib", "Rucaparib"),
  type   = c("APO", "Type II", "Type II",
             "Type III", "Type III", "Type III", "Type III"),
  stringsAsFactors = FALSE
)

load_rmsf <- function(sys_name) {
  f <- file.path(data_dir, paste0(sys_name, "_rmsf_residue.dat"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, skip = 1, col.names = c("resid", "rmsf"))
  df$system <- sys_name
  df
}

rmsf_list <- lapply(TYPE_MAP$system, load_rmsf)
names(rmsf_list) <- TYPE_MAP$system
rmsf_all <- bind_rows(rmsf_list[!sapply(rmsf_list, is.null)])
rmsf_all <- left_join(rmsf_all, TYPE_MAP, by = "system")

# Filter to CAT pocket residues
rmsf_pocket <- rmsf_all %>% filter(resid %in% CAT_POCKET)

# ---- Per-inhibitor pocket RMSF profile -------------------------------------
rmsf_pocket$label <- factor(rmsf_pocket$label,
                            levels = c("APO", "Talazoparib", "AZD5305",
                                       "Veliparib", "Niraparib", "Olaparib", "Rucaparib"))

p_profile <- ggplot(rmsf_pocket, aes(x = resid, y = rmsf, color = type, group = label)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(data = subset(rmsf_pocket, resid %in% CAT_FEATURES$resid),
             aes(shape = label), size = 1.5) +
  scale_color_manual(values = c("APO" = "#666666", "Type II" = "#2166AC", "Type III" = "#B2182B")) +
  scale_shape_manual(values = c(16, 17, 17, 15, 15, 15, 15)) +
  # Feature annotations
  annotate("rect", xmin = 859, xmax = 863, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3) +
  annotate("text", x = 861, y = max(rmsf_pocket$rmsf) * 0.95,
           label = "G-loop", size = 2, color = "grey40") +
  annotate("rect", xmin = 886, xmax = 890, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3) +
  annotate("text", x = 888, y = max(rmsf_pocket$rmsf) * 0.90,
           label = "D-loop", size = 2, color = "grey40") +
  annotate("rect", xmin = 903, xmax = 908, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3) +
  annotate("text", x = 905.5, y = max(rmsf_pocket$rmsf) * 0.95,
           label = "DON", size = 2, color = "grey40") +
  labs(x = "ART Residue", y = "RMSF (Å)", color = NULL,
       title = "CAT Pocket Residue Dynamics: Type II vs III") +
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        panel.grid.minor = element_blank())

# ---- Delta RMSF: Type II mean - Type III mean in pocket --------------------
pocket_delta <- rmsf_pocket %>%
  filter(type != "APO") %>%
  group_by(resid, type) %>%
  summarise(mean_rmsf = mean(rmsf), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = mean_rmsf) %>%
  mutate(
    delta = `Type II` - `Type III`,
    sign  = ifelse(delta > 0, "+", "-")
  ) %>%
  arrange(desc(abs(delta)))

cat("\n=== CAT Pocket Delta-RMSF (Type II - Type III) ===\n")
print(pocket_delta, n = 50)

# Key pocket residues with large differences
pocket_key <- pocket_delta %>% filter(abs(delta) > 0.3)

p_delta <- ggplot(pocket_key, aes(x = reorder(factor(resid), delta), y = delta, fill = sign)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("+" = "#B2182B", "-" = "#2166AC"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  labs(x = NULL, y = expression(Delta*"RMSF (Type II - Type III) [Å]"),
       title = "Key CAT Pocket Residues with Differential Dynamics") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank())

# ---- Export for PyMOL mapping ----
pocket_export <- pocket_delta %>%
  mutate(
    pymol_color = case_when(
      delta > 0.3 ~ "red",      # Type II more flexible
      delta < -0.3 ~ "blue",    # Type III more flexible  
      TRUE ~ "white"
    )
  )

# Generate PyMOL selection commands
cat("\n=== PyMOL commands for pocket mapping ===\n")
for (i in seq_len(nrow(pocket_export))) {
  if (pocket_export$pymol_color[i] != "white") {
    cat(sprintf("color %s, resi %d\n", pocket_export$pymol_color[i], pocket_export$resid[i]))
  }
}

write.csv(pocket_export, file.path(out_dir, "Fig5A_pocket_delta_rmsf.csv"), row.names = FALSE)

# ---- Final figure ----
library(patchwork)
master <- p_profile / p_delta +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "CAT Pocket Differential Dynamics: Type II vs III PARP1 Inhibitors",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

cairo_pdf(file.path(out_dir, "Fig5_Pocket_Dynamics.pdf"),
          width = 190/25.4, height = 140/25.4, pointsize = 7)
print(master)
dev.off()

cat("\nSaved: Fig5_Pocket_Dynamics.pdf\n")
cat("Saved: Fig5A_pocket_delta_rmsf.csv\n")

# ---- Summary ----
cat("\n========== POCKET MAPPING SUMMARY ==========\n")
cat(sprintf("Total CAT pocket residues analyzed: %d\n", length(CAT_POCKET)))
cat(sprintf("Residues with |ΔRMSF| > 0.3 Å: %d\n", nrow(pocket_key)))

# Which structural features show largest differences
pocket_key_feat <- left_join(pocket_key, CAT_FEATURES, by = "resid")
cat("\nKey differentially dynamic features:\n")
for (i in seq_len(min(10, nrow(pocket_key_feat)))) {
  r <- pocket_key_feat[i, ]
  feat <- ifelse(is.na(r$feature), "pocket", r$feature)
  cat(sprintf("  %d (%s): Δ=%.3f Å, feature=%s\n",
              r$resid, r$label %||% paste0("R", r$resid), r$delta, feat))
}

message("===== Pocket mapping complete =====")
