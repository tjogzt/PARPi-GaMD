#!/usr/bin/env Rscript
# 27-hbond_protein_dna.R — Protein-DNA H-bond occupancy analysis
#
# Uses cpptraj to compute H-bonds between PARP1 and DNA
# across all 7 S2 systems, then aggregates and visualizes.
#
# Residue layout (Amber 1-indexed):
#   Protein: :1-1015,1040-1041,1066,1070-1098
#   DNA:     :1016-1039,1042-1065
#   Zn:      :1067-1069
#   Water:   :1099-
#
# Output:
#   results/hbond_protein_dna/
#     hb_occupancy.csv       — per-residue-pair occupancy
#     hb_occupancy_heatmap.pdf — heatmap
#     hb_summary.csv         — per-system summary stats
library(tidyverse)
library(patchwork)

# ---- Config ----
DATA    <- "/Volumes/tjogzt4T/PARPi_data"
OUT_DIR <- file.path("results", "hbond_protein_dna")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CPPTRAJ <- "/opt/anaconda3/bin/cpptraj"

PRMTOP  <- file.path(DATA, "sys2_APO", "sys2_APO.prmtop")

systems <- c("sys2_APO", "sys2_AZD5305", "sys2_talazoparib",
             "sys2_olaparib", "sys2_niraparib", "sys2_rucaparib",
             "sys2_veliparib")

labels <- c("APO", "AZD5305", "Talazoparib", "Olaparib",
            "Niraparib", "Rucaparib", "Veliparib")

# ---- Step 1: Run cpptraj for each system ----
cat("\n=== Step 1: cpptraj H-bond detection ===\n")

for (i in seq_along(systems)) {
  sys <- systems[i]
  
  # Determine prmtop and DCD per system
  if (sys == "sys2_APO") {
    prmtop_file <- file.path(DATA, "sys2_APO", "sys2_APO.prmtop")
    dcd_file    <- file.path(DATA, "sys2_APO", "output.dcd")
    strip_extra <- "strip :WAT,Na+,Cl-"
  } else if (sys == "sys2_talazoparib") {
    prmtop_file <- file.path(DATA, "sys2_talazoparib", "sys2_talazoparib.prmtop")
    dcd_file    <- file.path(DATA, "sys2_talazoparib", "output.dcd")
    strip_extra <- "strip :WAT,Na+,Cl-,UNL"
  } else {
    prmtop_file <- PRMTOP
    dcd_file    <- file.path(DATA, sys, "output_stripped.dcd")
    strip_extra <- "strip :WAT,Na+,Cl-"
  }
  
  if (!file.exists(dcd_file)) {
    cat(sprintf("[SKIP] %s: no DCD\n", sys))
    next
  }
  
  sys_dir <- file.path(OUT_DIR, sys)
  dir.create(sys_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Skip if already done
  avg_file <- file.path(sys_dir, "hb_avg.dat")
  if (file.exists(avg_file) && file.info(avg_file)$size > 10) {
    cat(sprintf("  %s... (cached)\n", sys))
    next
  }
  
  cpptraj_in <- file.path(sys_dir, "hb.in")
  hbond_cmd <- "hbond HB :1-1015,1040-1041,1066,1070-1098 acceptormask :1016-1039,1042-1065 out hb.dat avgout hb_avg.dat dist 3.5 angle 120.0 nointramol"
  writeLines(c(
    sprintf("parm %s", normalizePath(prmtop_file)),
    sprintf("trajin %s", normalizePath(dcd_file)),
    strip_extra,
    hbond_cmd,
    "run",
    "quit"
  ), cpptraj_in)
  
  cat(sprintf("  %s...", sys))
  cmd <- sprintf("cd %s && %s -i hb.in", sys_dir, CPPTRAJ)
  ret <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
  cat(" done\n")
}

cat("cpptraj complete.\n")

# ---- Step 2: Parse per-system H-bond averages ----
cat("\n=== Step 2: Parse H-bond data ===\n")

# hb_avg.dat format (space-delimited):
# #Acceptor              DonorH          Donor   Frames         Frac      AvgDist       AvgAng
# DG_1039@OP1         SER_16@HG      SER_16@OG      194       0.9949       2.6962     165.1895

all_hb <- list()

for (i in seq_along(systems)) {
  sys <- systems[i]
  avg_file <- file.path(OUT_DIR, sys, "hb_avg.dat")
  
  if (!file.exists(avg_file)) next
  
  lines <- readLines(avg_file)
  if (length(lines) < 2) next  # only header or empty
  
  # Read as space-delimited table; skip header line (#Acceptor ...)
  lines <- readLines(avg_file)
  if (length(lines) <= 1) next  # only header or empty
  
  # Remove leading # from header and parse
  header_line <- sub("^#", "", lines[1])
  header_fields <- strsplit(trimws(header_line), "[[:space:]]+")[[1]]
  
  dat <- read.table(text = lines[-1], header = FALSE,
                    col.names = header_fields,
                    stringsAsFactors = FALSE)
  
  if (nrow(dat) == 0) next
  
  # Parse residue info from column names
  # Acceptor: DG_1039@OP1 → dna_resname=DG, dna_resnum=1039, dna_atom=OP1
  parse_label <- function(x) {
    m <- regmatches(x, regexec("^([A-Z0-9]+)_([0-9]+)@(.+)$", x))[[1]]
    if (length(m) == 4) {
      data.frame(resname = m[2], resnum = as.integer(m[3]), atom = m[4],
                 stringsAsFactors = FALSE)
    } else {
      data.frame(resname = NA_character_, resnum = NA_integer_, 
                 atom = NA_character_, stringsAsFactors = FALSE)
    }
  }
  
  acceptor <- do.call(rbind, lapply(dat$Acceptor, parse_label))
  names(acceptor) <- c("dna_resname", "dna_resnum", "dna_atom")
  
  donor <- do.call(rbind, lapply(dat$Donor, parse_label))
  names(donor) <- c("prot_resname", "prot_resnum", "prot_atom")
  
  hb_pairs <- data.frame(
    system       = sys,
    dna_resname  = acceptor$dna_resname,
    dna_resnum   = acceptor$dna_resnum,
    dna_atom     = acceptor$dna_atom,
    prot_resname = donor$prot_resname,
    prot_resnum  = donor$prot_resnum,
    prot_atom    = donor$prot_atom,
    occupancy    = dat$Frac,
    frames       = dat$Frames,
    avg_distance = dat$AvgDist,
    avg_angle    = dat$AvgAng,
    stringsAsFactors = FALSE
  )
  
  if (nrow(hb_pairs) > 0) {
    all_hb[[sys]] <- hb_pairs
  }
}

hb_df <- bind_rows(all_hb)
hb_df$system <- factor(hb_df$system, levels = systems, labels = labels)

cat(sprintf("Total H-bond pairs detected: %d\n", nrow(hb_df)))

# ---- Step 3: Per-system summary ----
cat("\n=== Step 3: Summary statistics ===\n")

summary_df <- hb_df %>%
  group_by(system) %>%
  summarise(
    n_hbonds        = n(),
    n_prot_residues = n_distinct(prot_resnum),
    n_dna_residues  = n_distinct(dna_resnum),
    mean_occupancy  = mean(occupancy),
    sd_occupancy    = sd(occupancy),
    n_strong_bonds  = sum(occupancy > 0.5),
    n_weak_bonds    = sum(occupancy > 0.1 & occupancy <= 0.5),
    .groups = "drop"
  ) %>%
  arrange(desc(n_strong_bonds))

print(summary_df)

write.csv(summary_df, file.path(OUT_DIR, "hb_summary.csv"), row.names = FALSE)

# ---- Step 4: Occupancy heatmap ----
cat("\n=== Step 4: Heatmap ===\n")

# Pivot to matrix: rows=protein residues, cols=DNA residues
# Take max occupancy per residue pair (multiple atom pairs may exist)
hb_res <- hb_df %>%
  group_by(system, prot_resnum, dna_resnum) %>%
  summarise(occupancy = max(occupancy), .groups = "drop")

# Compute average occupancy across all systems for ordering
avg_occ <- hb_res %>%
  group_by(prot_resnum, dna_resnum) %>%
  summarise(avg_occ = mean(occupancy), .groups = "drop")

# Order protein residues by total DNA H-bond occupancy
prot_order <- avg_occ %>%
  group_by(prot_resnum) %>%
  summarise(total = sum(avg_occ), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(prot_resnum)

dna_order <- avg_occ %>%
  group_by(dna_resnum) %>%
  summarise(total = sum(avg_occ), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(dna_resnum)

hb_res$prot_resnum <- factor(hb_res$prot_resnum, levels = prot_order)
hb_res$dna_resnum  <- factor(hb_res$dna_resnum, levels = dna_order)

# Per-system heatmaps
plot_list <- list()
for (sys_name in labels) {
  sub <- hb_res %>% filter(system == sys_name)
  p <- ggplot(sub, aes(x = dna_resnum, y = prot_resnum, fill = occupancy)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                         name = "Occupancy") +
    labs(title = sys_name, x = "DNA Residue", y = "Protein Residue") +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = element_text(size = 5),
      plot.title = element_text(size = 8, face = "bold"),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm")
    )
  plot_list[[sys_name]] <- p
}

# Combine: 2 columns x 4 rows
combined <- wrap_plots(plot_list, ncol = 2)

ggsave(file.path(OUT_DIR, "hb_occupancy_heatmap.pdf"),
       combined, width = 8, height = 14, device = cairo_pdf)

# Also a combined average heatmap
p_avg <- ggplot(avg_occ %>%
         mutate(prot_resnum = factor(prot_resnum, levels = prot_order),
                dna_resnum  = factor(dna_resnum, levels = dna_order)),
       aes(x = dna_resnum, y = prot_resnum, fill = avg_occ)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 0.8),
                       name = "Mean\nOccupancy") +
  labs(title = "Average Protein-DNA H-bond Occupancy",
       x = "DNA Residue", y = "Protein Residue") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(file.path(OUT_DIR, "hb_occupancy_avg_heatmap.pdf"),
       p_avg, width = 7, height = 8, device = cairo_pdf)

# ---- Step 5: Bar chart of strong H-bonds per system ----
cat("\n=== Step 5: Bar chart ===\n")

p_bar <- ggplot(summary_df, aes(x = reorder(system, n_strong_bonds), 
                                 y = n_strong_bonds)) +
  geom_col(aes(fill = n_strong_bonds), width = 0.6) +
  geom_text(aes(label = n_strong_bonds), hjust = -0.3, size = 2.5) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = NULL, y = "Strong H-bonds (occupancy > 50%)",
       title = "Protein-DNA H-bond Network Size") +
  coord_flip() +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "hb_strong_bonds_bar.pdf"),
       p_bar, width = 5, height = 3.5, device = cairo_pdf)

# ---- Step 6: Top residue pairs per system ----
cat("\n=== Step 6: Top H-bond pairs ===\n")

top_pairs <- hb_df %>%
  filter(occupancy > 0.5) %>%
  group_by(system) %>%
  slice_max(order_by = occupancy, n = 10) %>%
  arrange(system, desc(occupancy)) %>%
  select(system, prot_resname, prot_resnum, prot_atom,
         dna_resname, dna_resnum, dna_atom, occupancy)

print(top_pairs, n = 100)
write.csv(top_pairs, file.path(OUT_DIR, "hb_top_pairs.csv"), row.names = FALSE)

# ---- Step 7: Occupancy CSV ----
occ_csv <- hb_res %>%
  pivot_wider(names_from = system, values_from = occupancy, 
              values_fill = 0)
write.csv(occ_csv, file.path(OUT_DIR, "hb_occupancy.csv"), row.names = FALSE)

cat(sprintf("\n✅ All outputs in %s/\n", OUT_DIR))
cat("  hb_summary.csv, hb_top_pairs.csv, hb_occupancy.csv\n")
cat("  hb_occupancy_heatmap.pdf, hb_occupancy_avg_heatmap.pdf\n")
cat("  hb_strong_bonds_bar.pdf\n")
