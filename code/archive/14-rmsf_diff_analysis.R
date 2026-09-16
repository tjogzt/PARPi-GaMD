#!/usr/bin/env Rscript
# 14-rmsf_diff_analysis.R — Per-domain RMSF comparison: APO vs talazoparib (Type II)
# Identifies residues with significant delta-RMSF to map allosteric pathway
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---- Config ----------------------------------------------------------------
data_dir  <- "results/analysis"
out_dir   <- "results/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# PARP1 domain boundaries (based on PDB 4DQY / UniProt P09874)
domain_map <- data.frame(
  domain   = c("ZnF1", "ZnF2", "ZnF3", "BRCT", "WGR", "HD", "ART"),
  start    = c(1, 97, 215, 384, 518, 662, 788),
  end      = c(96, 214, 383, 517, 661, 787, 1011),
  color    = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"),
  stringsAsFactors = FALSE
)

# ---- Theme -----------------------------------------------------------------
theme_7pt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"), strip.background = element_rect(fill = "grey95"))

# ---- Load RMSF data --------------------------------------------------------
read_rmsf <- function(sys, dom) {
  f <- file.path(data_dir, paste0(sys, "_rmsf_", dom, ".dat"))
  if (!file.exists(f)) return(NULL)
  df <- read.table(f, skip = 1, col.names = c("resid", "rmsf"))
  df$system <- ifelse(grepl("APO", sys), "APO", "Talazoparib")
  df$domain <- dom
  df
}

sys_list <- c("sys2_APO", "sys2_talazoparib")
dom_list <- c("all", "HD", "ART", "DNA")

rmsf_all <- do.call(rbind, unlist(lapply(sys_list, function(s) {
  lapply(dom_list, function(d) read_rmsf(s, d))
}), recursive = FALSE))

rmsf_all$system <- factor(rmsf_all$system, levels = c("APO", "Talazoparib"))
rmsf_all$domain <- factor(rmsf_all$domain, levels = c("all", "HD", "ART", "DNA"))

# ---- Panel A: Per-domain RMSF overlay --------------------------------------
p_a <- ggplot(rmsf_all, aes(x = resid, y = rmsf, color = system)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~ domain, scales = "free", ncol = 2) +
  scale_color_manual(values = c("APO" = "#2166AC", "Talazoparib" = "#B2182B")) +
  labs(x = "Residue", y = "RMSF (Å)", color = NULL,
       title = "A  Per-Domain RMSF: APO vs Talazoparib") +
  theme_7pt + theme(legend.position = "bottom")

# ---- Panel B: Delta-RMSF (Talazoparib - APO) for HD domain -----------------
rmsf_hd <- rmsf_all %>% filter(domain == "HD") %>%
  select(resid, rmsf, system) %>%
  pivot_wider(names_from = system, values_from = rmsf) %>%
  mutate(delta = Talazoparib - APO,
         sign  = ifelse(delta > 0, "+", "-"),
         abs_delta = abs(delta))

# Top 20 most changed residues
top20 <- rmsf_hd %>% arrange(desc(abs_delta)) %>% head(20)
cat("\n=== Top 20 HD residues with largest |ΔRMSF| (Talazoparib - APO) ===\n")
print(top20, n = 20)

# Panel B: Delta RMSF bar plot (HD domain only)
p_b <- ggplot(rmsf_hd, aes(x = resid, y = delta, fill = sign)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("+" = "#B2182B", "-" = "#2166AC"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  # Highlight CAT pocket proximal region (HD residues near active site ~662-700)
  annotate("rect", xmin = 662, xmax = 710, ymin = -Inf, ymax = Inf, 
           fill = "yellow", alpha = 0.12) +
  annotate("text", x = 686, y = max(rmsf_hd$delta, na.rm=TRUE)*0.9,
           label = "CAT\npocket\nproximal", size = 2.2, color = "grey40", hjust = 0.5) +
  labs(x = "HD Domain Residue", y = expression(Delta*"RMSF (Talazoparib - APO) [Å]"),
       title = "B  HD Domain: Talazoparib-Induced Rigidification") +
  theme_7pt

# ---- Panel C: RMSF correlation with S1 well_depth --------------------------
# Summary metrics per domain
rmsf_summary <- rmsf_all %>%
  group_by(system, domain) %>%
  summarise(
    mean_rmsf  = mean(rmsf, na.rm = TRUE),
    sd_rmsf    = sd(rmsf, na.rm = TRUE),
    max_rmsf   = max(rmsf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = system, values_from = c(mean_rmsf, sd_rmsf, max_rmsf))

cat("\n=== RMSF Summary by Domain ===\n")
print(as.data.frame(rmsf_summary))

# ---- Panel C: HD domain RMSF heatmap-style (binned) ------------------------
# Bin residues into 20-residue windows
rmsf_hd_binned <- rmsf_hd %>%
  mutate(bin = floor((resid - 662) / 15) * 15 + 662) %>%
  group_by(bin) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = paste0(bin, "-", bin + 14))

p_c <- ggplot(rmsf_hd_binned, aes(x = reorder(label, bin), y = 1, fill = mean_delta)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       name = expression(Delta*"RMSF (Å)"), midpoint = 0) +
  labs(x = "HD Domain Residue Range", y = NULL,
       title = "C  Binned ΔRMSF: Talazoparib Rigidifies HD Domain") +
  theme_7pt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

# ---- Assemble figure -------------------------------------------------------
master <- (p_a / (p_b | p_c)) +
  plot_layout(heights = c(2, 1.2)) +
  plot_annotation(
    title = "PARP1 Trapping: Talazoparib Induces HD Domain Rigidification",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

cairo_pdf(file.path(out_dir, "Fig_RMSF_Analysis.pdf"), width = 190/25.4, height = 160/25.4, pointsize = 7)
print(master)
dev.off()

cat("\nSaved: Fig_RMSF_Analysis.pdf\n")

# ---- Key finding summary ---------------------------------------------------
cat("\n========== KEY RMSF FINDINGS ==========\n")
# Mean delta per domain
for (dom in c("HD", "ART", "DNA", "all")) {
  dd <- rmsf_all %>% filter(domain == dom) %>%
    select(resid, rmsf, system) %>%
    pivot_wider(names_from = system, values_from = rmsf) %>%
    mutate(delta = Talazoparib - APO)
  cat(sprintf("%s: mean ΔRMSF = %.3f Å, %s\n", dom, mean(dd$delta, na.rm=TRUE),
              ifelse(mean(dd$delta) < 0, "(more rigid with talazoparib)", "(more flexible)")))
}

# Count significantly rigidified vs flexibilized residues in HD
hd_dd <- rmsf_all %>% filter(domain == "HD") %>%
  select(resid, rmsf, system) %>%
  pivot_wider(names_from = system, values_from = rmsf) %>%
  mutate(delta = Talazoparib - APO)

n_rigid <- sum(hd_dd$delta < -0.5, na.rm = TRUE)
n_flex  <- sum(hd_dd$delta > 0.5, na.rm = TRUE)
cat(sprintf("HD domain: %d residues rigidified (Δ<-0.5), %d flexibilized (Δ>0.5)\n", n_rigid, n_flex))

message("===== RMSF analysis complete =====")
