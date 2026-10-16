# 11b_s5_composite_figure.R — Figure S5 composite (S1/S2 HD-ART + AAI)
#
# Purpose:  Render the Figure S5 composite: seven S1-vs-S2 HD-ART facet panels
#           in a 4+3 square grid (S1/S2 legend inset into the vacant eighth
#           cell) with (B) the S1/S2 span ratio and (C) the S1 well-depth map
#           as square panels below.
# Created:  2026-10-14
# Depends:  dplyr, ggplot2, ggrepel, patchwork
# Inputs:   data/analysis/sys1_*_pmf_c3.xvg ; data/analysis/pmf-c3-sys2_*_CV2_cv.dat.xvg ;
#           data/analysis/r6/r6_panel_values.csv ; common/ligands.csv
# Run:      Rscript code/11b_s5_composite_figure.R   (from the repository root)
args_h <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args_h[grep("^--file=", args_h)])))
source(file.path(script_dir, "..", "common", "helpers.R"))
source(file.path(script_dir, "..", "common", "ligands.R"))
library(ggplot2); library(dplyr); library(patchwork)

data_dir <- "data/analysis"
out_dir  <- "figures/pdf"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sys_meta <- data.frame(
  ligand   = c("APO", "AZD5305", "olaparib", "talazoparib", "veliparib", "niraparib", "rucaparib"),
  label    = c("APO", "AZD5305", "Olaparib", "Talazoparib", "Veliparib", "Niraparib", "Rucaparib"),
  class    = c("APO", "Unknown", "Type_II", "Type_II", "Type_III", "Type_III", "Type_III"),
  color    = c("grey40", "darkorange", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3"),
  stringsAsFactors = FALSE
)

# ---- Load S1 (CAT-only) and S2 (CV2) HD-ART PMFs ----------------------------
s1_pmf <- list(); s2_pmf <- list()
for (i in seq_len(nrow(sys_meta))) {
  lig <- sys_meta$ligand[i]
  f1 <- file.path(data_dir, sprintf("sys1_%s_pmf_c3.xvg", lig))
  if (file.exists(f1)) {
    df <- read_pmf(f1); df$ligand <- lig; df$label <- sys_meta$label[i]
    df$system <- "S1: CAT-only"; s1_pmf[[lig]] <- df
  }
  f2 <- file.path(data_dir, sprintf("pmf-c3-sys2_%s_CV2_cv.dat.xvg", lig))
  if (file.exists(f2)) {
    df <- read_pmf(f2); df$ligand <- lig; df$label <- sys_meta$label[i]
    df$system <- "S2: DNA-bound"; s2_pmf[[lig]] <- df
  }
}
pmf_all <- bind_rows(bind_rows(s1_pmf), bind_rows(s2_pmf))
pmf_all$system <- factor(pmf_all$system, levels = c("S1: CAT-only", "S2: DNA-bound"))

theme_pmf <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    plot.margin       = ggplot2::margin(1, 1, 1, 3, unit = "mm"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 9, face = "bold"),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 8),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 8),
    legend.key.size   = unit(0.3, "cm"),
    strip.text        = element_text(size = 8, face = "bold"),
    strip.background  = element_rect(fill = "grey95")
  )

# ---- Facet: seven square panels (4+3) with legend inset ---------------------
p_facet <- ggplot(pmf_all, aes(x = RC, y = PMF_norm, color = system)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("S1: CAT-only" = "#2166AC",
                                "S2: DNA-bound" = "#B2182B")) +
  facet_wrap(~ label, ncol = 4, scales = "free_x") +
  labs(x = "HD-ART COM Distance (Å)", y = "Free Energy (kcal/mol)", color = NULL) +
  theme_pmf +
  theme(legend.position = "none", aspect.ratio = 1)

p_facet_leg <- p_facet + theme(legend.position = "bottom", aspect.ratio = NULL,
                               legend.key.width = unit(0.5, "cm")) +
  guides(color = guide_legend(ncol = 1))
leg_grob <- {
  g <- ggplot2::ggplotGrob(p_facet_leg)
  g$grobs[[which(vapply(g$grobs, function(x) x$name, "") == "guide-box")]]
}

# ---- (B) S1/S2 span ratio (code 12 verbatim core) ---------------------------
pv2 <- read.csv(file.path(data_dir, "r6", "r6_panel_values.csv"))
base <- pv2[!grepl("_prot", pv2$ligand), ]
prot <- pv2[grepl("_prot", pv2$ligand), ]
sc <- rbind(
  data.frame(x = base$S2CV2_C3, y = base$S1_C3, label = base$label, filled = TRUE),
  data.frame(x = prot$S2CV2_C3, y = prot$S1_C3, label = prot$label, filled = FALSE))
sc$fill_key <- ifelse(sc$filled, "base", "prot")
p_ratio <- ggplot(sc, aes(x = x, y = y, label = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = label, shape = fill_key), size = 3.0, stroke = 0.9) +
  scale_shape_manual(values = c("base" = 16, "prot" = 1), guide = "none") +
  ggrepel::geom_text_repel(size = 3.0, force = 2, box.padding = 0.35, max.overlaps = Inf,
                           seed = 49, family = "Arial") +
  scale_color_manual(values = c("APO" = "grey40", "AZD5305" = "darkorange",
                                "Olaparib" = "#E41A1C", "Talazoparib" = "#FF7F00",
                                "Veliparib" = "#377EB8", "Niraparib" = "#4DAF4A",
                                "Rucaparib" = "#984EA3", "Fluzoparib" = "#C23531",
                                "Pamiparib" = "#3D6BA8", "Senaparib" = "#9D2933",
                                "Niraparib (prot.)" = "#4DAF4A",
                                "Rucaparib (prot.)" = "#984EA3",
                                "Veliparib (prot.)" = "#377EB8")) +
  labs(x = "S2 CV2 span (kcal/mol, C3)",
       y = "S1 span (kcal/mol, C3)") +
  theme_7pt + theme(legend.position = "none",
                    plot.margin = ggplot2::margin(1, 1, 1, 3, unit = "mm"),
                    aspect.ratio = 1)

# ---- (C) 2D well-depth map (code 13 verbatim core) --------------------------
sys10 <- load_ligands()
sys10$type <- gsub("_", " ", sys10$class)
pv <- read.csv(file.path(data_dir, "r6", "r6_panel_values.csv"))
b8 <- c("azd5305", "olaparib", "talazoparib", "veliparib", "rucaparib",
        "fluzoparib", "pamiparib", "senaparib")
m8 <- match(b8, pv$ligand)
m10 <- match(b8, tolower(sys10$ligand))   # ligands.csv uses "AZD5305" (upper case)
plot_df2 <- data.frame(
  ligand = b8, label = sys10$label[m10], type = sys10$type[m10],
  wd_ratio = pv$AAI_C3[m8], wd_S1 = pv$S1_C3[m8],
  state = "base", shape2 = sys10$type[m10], stringsAsFactors = FALSE)
nira_neut <- data.frame(ligand = "niraparib", label = "Niraparib\n(neutral)",
                        type = "Type III", wd_ratio = pv$AAI_C3[pv$ligand == "niraparib"],
                        wd_S1 = pv$S1_C3[pv$ligand == "niraparib"],
                        state = "nira_neutral", shape2 = "nira", stringsAsFactors = FALSE)
nira_prot <- data.frame(ligand = "niraparib", label = "Niraparib\n(protonated)",
                        type = "Type III", wd_ratio = pv$AAI_C3[pv$ligand == "niraparib_prot"],
                        wd_S1 = pv$S1_C3[pv$ligand == "niraparib_prot"],
                        state = "nira_protonated", shape2 = "nira", stringsAsFactors = FALSE)
ruca_prot2 <- data.frame(ligand = "rucaparib", label = "Rucaparib\n(protonated)",
                         type = "Type III", wd_ratio = pv$AAI_C3[pv$ligand == "rucaparib_prot"],
                         wd_S1 = pv$S1_C3[pv$ligand == "rucaparib_prot"],
                         state = "ruca_protonated", shape2 = "ruca", stringsAsFactors = FALSE)
# (2026-10-15) veliparib protonated state added: the audited-dominant panel
# comprises all three state-corrected ligands; this point was missing.
veli_prot2 <- data.frame(ligand = "veliparib", label = "Veliparib\n(protonated)",
                         type = "Type III", wd_ratio = pv$AAI_C3[pv$ligand == "veliparib_prot"],
                         wd_S1 = pv$S1_C3[pv$ligand == "veliparib_prot"],
                         state = "veli_protonated", shape2 = "veli", stringsAsFactors = FALSE)
plot_df2 <- bind_rows(plot_df2, nira_neut, nira_prot, ruca_prot2, veli_prot2)
p_2d <- ggplot(plot_df2, aes(x = wd_ratio, y = wd_S1)) +
  geom_point(aes(color = type, shape = shape2), size = 4) +
  scale_shape_manual(values = c("Type II" = 16, "Type III" = 17,
                                "Extension" = 18, "Unknown" = 15,
                                "nira" = 1, "ruca" = 1, "veli" = 1),
                     breaks = c("Type II", "Type III", "Extension", "Unknown"),
                     name = "Classification") +
  ggrepel::geom_text_repel(aes(label = label), size = 3.0, max.overlaps = Inf,
                           seed = 49, min.segment.length = 0.3, box.padding = 0.6,
                           force = 4, max.iter = 5000, family = "Arial") +
  scale_color_manual(values = c("Type II" = "#E41A1C", "Type III" = "#377EB8",
                                "Extension" = "#C23531", "Unknown" = "darkorange"),
                     breaks = c("Type II", "Type III", "Extension", "Unknown"),
                     name = "Classification") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "S1/S2 span ratio",
       y = "S1 HD--ART span (kcal/mol)") +
  theme_7pt +
  theme(legend.position = "none",
        plot.margin = ggplot2::margin(1, 4, 1, 3, unit = "mm"),
        aspect.ratio = 1)

# ---- Combined legend cell (row 2, column 4): S1/S2 + Classification --------
# Both legends are gathered into the vacant eighth cell of the facet grid so
# that (C) renders as a plain square panel and the two panel rows keep the
# same total width as the (B)|(C) row.
p_2d_leg <- p_2d + theme(legend.position = "bottom",
                         legend.direction = "horizontal",
                         legend.spacing.x = unit(0.15, "cm")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
         shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))
leg_class <- {
  g <- ggplot2::ggplotGrob(p_2d_leg)
  g$grobs[[which(vapply(g$grobs, function(x) x$name, "") == "guide-box")]]
}
leg_cell <- wrap_plots(leg_grob, leg_class, ncol = 1)
p_facet_final <- p_facet + inset_element(leg_cell,
                                         left = 0.762, bottom = -0.122,
                                         right = 0.992, top = 0.308)

# ---- Compose: facets on top, (B)|(C) row below ------------------------------
row3 <- p_ratio | p_2d
final <- p_facet_final / row3 + plot_layout(heights = c(1.3, 1))

cairo_pdf(file.path(out_dir, "Fig_S1S2_spanratio_composite.pdf"),
          width = 178/25.4, height = 195/25.4, pointsize = 8)
print(final)
dev.off()
message("Saved: Fig_S1S2_spanratio_composite.pdf")
