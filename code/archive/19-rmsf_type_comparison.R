#!/usr/bin/env Rscript
# 19-rmsf_type_comparison.R — Per-residue RMSF: Type II vs Type III comparison
# Generates Fig 3D (RMSF diff) and Fig 5B (CAT pocket mapping)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# PARP1 domain boundaries (UniProt P09874, S2 prmtop residue numbering)
DOMAIN_MAP <- data.frame(
  domain = c("ZnF1", "ZnF2", "ZnF3", "BRCT", "WGR", "HD", "ART"),
  start  = c(1, 97, 215, 384, 518, 662, 788),
  end    = c(96, 214, 383, 517, 661, 787, 1014),
  color  = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"),
  stringsAsFactors = FALSE
)

# Inhibitor classification
TYPE_MAP <- data.frame(
  system = c("sys2_APO", "sys2_AZD5305", "sys2_talazoparib",
             "sys2_veliparib", "sys2_niraparib", "sys2_olaparib", "sys2_rucaparib"),
  label  = c("APO", "AZD5305", "Talazoparib",
             "Veliparib", "Niraparib", "Olaparib", "Rucaparib"),
  type   = c("APO", "Type II", "Type II",
             "Type III", "Type III", "Type III", "Type III"),
  color  = c("#666666", "#2166AC", "#2166AC",
             "#B2182B", "#B2182B", "#B2182B", "#B2182B"),
  stringsAsFactors = FALSE
)

# ---- Theme -----------------------------------------------------------------
theme_7pt <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_rect(fill = "grey95")
  )

# ---- Load RMSF data --------------------------------------------------------
load_rmsf <- function(sys_name) {
  f <- file.path(data_dir, paste0(sys_name, "_rmsf_residue.dat"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, skip = 1, col.names = c("resid", "rmsf"))
  df$system <- sys_name
  df
}

all_systems <- TYPE_MAP$system
rmsf_list <- lapply(all_systems, load_rmsf)
names(rmsf_list) <- all_systems
rmsf_list <- rmsf_list[!sapply(rmsf_list, is.null)]

if (length(rmsf_list) < 3) {
  stop("Insufficient RMSF data. Run 19-rmsf_all_systems.py first.")
}

rmsf_all <- bind_rows(rmsf_list)
rmsf_all <- left_join(rmsf_all, TYPE_MAP, by = "system")
rmsf_all$type <- factor(rmsf_all$type, levels = c("APO", "Type II", "Type III"))

cat(sprintf("Loaded RMSF data for %d systems\n", length(unique(rmsf_all$system))))
cat("Systems:", paste(unique(rmsf_all$label), collapse = ", "), "\n")

# ---- Domain annotation -----------------------------------------------------
annotate_domain <- function(resid) {
  for (i in seq_len(nrow(DOMAIN_MAP))) {
    if (resid >= DOMAIN_MAP$start[i] && resid <= DOMAIN_MAP$end[i]) {
      return(DOMAIN_MAP$domain[i])
    }
  }
  return(NA_character_)
}
rmsf_all$domain <- sapply(rmsf_all$resid, annotate_domain)
rmsf_all$domain <- factor(rmsf_all$domain, levels = DOMAIN_MAP$domain)

# ---- Panel A: RMSF overlay by domain (Type II vs III vs APO) --------------
# Focus on HD + ART domains
rmsf_catalytic <- rmsf_all %>% filter(domain %in% c("HD", "ART"))

p_a <- ggplot(rmsf_catalytic, aes(x = resid, y = rmsf, color = type, group = label)) +
  geom_line(linewidth = 0.25, alpha = 0.7) +
  facet_wrap(~ domain, scales = "free_x", ncol = 1, strip.position = "right") +
  scale_color_manual(values = c("APO" = "#666666", "Type II" = "#2166AC", "Type III" = "#B2182B")) +
  labs(x = "Residue", y = "RMSF (Å)", color = NULL,
       title = "A  Per-Residue RMSF: HD + ART Domains") +
  theme_7pt + theme(legend.position = "bottom")

# ---- Panel B: Delta-RMSF (Type II mean - Type III mean) -------------------
# Compute mean RMSF per residue per type
rmsf_type_mean <- rmsf_all %>%
  filter(domain %in% c("HD", "ART")) %>%
  group_by(resid, domain, type) %>%
  summarise(mean_rmsf = mean(rmsf, na.rm = TRUE), .groups = "drop") %>%
  filter(type != "APO") %>%
  pivot_wider(names_from = type, values_from = mean_rmsf) %>%
  mutate(
    delta = `Type II` - `Type III`,
    sign  = ifelse(delta > 0, "+", "-"),
    abs_delta = abs(delta)
  )

# Mark CAT pocket region (ART residues near active site: 859-908 based on PDB 4DQY)
# and HD proximal region (residues near ART interface)
cat_sig <- rmsf_type_mean %>%
  filter(domain == "ART") %>%
  arrange(desc(abs_delta)) %>%
  head(20)

cat(sprintf("\n=== Top 20 ART residues with largest |ΔRMSF| (Type II - Type III) ===\n"))
print(cat_sig, n = 20)

p_b <- ggplot(rmsf_type_mean, aes(x = resid, y = delta, fill = sign)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("+" = "#B2182B", "-" = "#2166AC"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  facet_wrap(~ domain, scales = "free_x", ncol = 1) +
  # Highlight CAT pocket
  annotate("rect", xmin = 859, xmax = 908, ymin = -Inf, ymax = Inf,
           fill = "yellow", alpha = 0.08) +
  labs(x = "Residue", y = expression(Delta*"RMSF (Type II - Type III) [Å]"),
       title = "B  Differential RMSF: Type II vs Type III") +
  theme_7pt

# ---- Panel C: CAT pocket zoom (ART residues 850-920) -----------------------
rmsf_cat <- rmsf_type_mean %>% filter(resid >= 850, resid <= 920)

# Structural annotations for key ART residues (based on 4DQY)
cat_features <- data.frame(
  resid = c(862, 863, 872, 877, 878, 888, 890, 896, 903, 904, 907, 908),
  label = c("G862", "R863", "S872", "G877", "R878", "S888", "D890", "Y896", "Y903", "S904", "Y907", "E908"),
  feature = c("loop", "HB", "HB", "loop", "HB", "loop", "HB", "HB", "loop", "loop", "loop", "HB"),
  stringsAsFactors = FALSE
)

# Simple label for key residues
cat_labels <- merge(cat_features, rmsf_cat, by = "resid")
cat_labels <- cat_labels[abs(cat_labels$delta) > 0.1, ]

p_c <- ggplot(rmsf_cat, aes(x = resid, y = delta, fill = sign)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("+" = "#B2182B", "-" = "#2166AC"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_text(data = cat_labels,
            aes(x = resid, y = delta, label = label),
            size = 1.8, vjust = -1, angle = 45, hjust = 0, inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(850, 920, 10)) +
  labs(x = "ART Domain Residue", y = expression(Delta*"RMSF (Å)"),
       title = "C  CAT Pocket Region: Type II Rigidification") +
  theme_7pt + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

# ---- Panel D: Per-inhibitor HD vs ART scatter ------------------------------
rmsf_dom_summary <- rmsf_all %>%
  filter(domain %in% c("HD", "ART")) %>%
  group_by(system, label, type, domain) %>%
  summarise(
    mean_rmsf = mean(rmsf, na.rm = TRUE),
    sd_rmsf   = sd(rmsf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = domain, values_from = c(mean_rmsf, sd_rmsf))

p_d <- ggplot(rmsf_dom_summary, aes(x = mean_rmsf_HD, y = mean_rmsf_ART, color = type)) +
  geom_point(size = 2.5) +
  geom_text(aes(label = label), size = 2, vjust = -1, hjust = 0.5) +
  scale_color_manual(values = c("APO" = "#666666", "Type II" = "#2166AC", "Type III" = "#B2182B")) +
  labs(x = "HD Mean RMSF (Å)", y = "ART Mean RMSF (Å)",
       title = "D  HD vs ART Rigidity Landscape") +
  theme_7pt + theme(legend.position = "none")

# ---- Assemble Figure -------------------------------------------------------
master <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1, 1), heights = c(1.5, 1)) +
  plot_annotation(
    title = "PARP1 Trapping: Type II Inhibitors Induce HD-ART Rigidification",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

cairo_pdf(file.path(out_dir, "Fig_RMSF_Type_Comparison.pdf"),
          width = 190/25.4, height = 180/25.4, pointsize = 7)
print(master)
dev.off()
cat("\nSaved: Fig_RMSF_Type_Comparison.pdf\n")

# ---- Statistical summary ---------------------------------------------------
cat("\n========== SUMMARY ==========\n")

# Per-domain mean RMSF by type
dom_summary <- rmsf_all %>%
  filter(domain %in% c("HD", "ART")) %>%
  group_by(domain, type) %>%
  summarise(
    mean = mean(rmsf, na.rm = TRUE),
    sd   = sd(rmsf, na.rm = TRUE),
    n    = n(),
    .groups = "drop"
  )
print(as.data.frame(dom_summary))

# Key finding: Type II vs III difference
hddiff <- rmsf_type_mean %>%
  group_by(domain) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    n_rigidified = sum(delta < -0.1, na.rm = TRUE),
    n_flexibilized = sum(delta > 0.1, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nType II vs Type III delta-RMSF by domain:\n")
print(as.data.frame(hddiff))

message("===== RMSF Type comparison complete =====")
