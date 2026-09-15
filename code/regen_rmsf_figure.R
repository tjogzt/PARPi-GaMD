#!/usr/bin/env Rscript
# Regenerate Fig_S2_RMSF_All.pdf — unified-protocol RMSF (CAT-domain alignment, production, 50 ps)
# (A) HD per-residue RMSF profiles  (B) domain-mean bars  (C) Type II-III dRMSF
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

theme_set(theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right", legend.key.size = unit(0.3, "cm")))

# China-style palette
cols <- c(APO = "#7A7A7A", talazoparib = "#C23531", olaparib = "#3D6BA8",
          niraparib = "#9D2933", rucaparib = "#177CB0",
          veliparib = "#B36B2E", AZD5305 = "#5E8C5E")

TYPE_MAP <- data.table(system = c("talazoparib", "olaparib", "niraparib", "rucaparib", "veliparib"),
                       type = c("II", "II", "III", "III", "III"))
HD <- 662:787; ART <- 788:1014

dir <- "results/analysis/rmsf_recomp"
sys_names <- c("APO", "talazoparib", "olaparib", "niraparib", "rucaparib", "veliparib", "AZD5305")

load_rmsf <- function(s) {
  f <- file.path(dir, paste0(s, "_rmsf_residue.dat"))
  x <- fread(f, skip = 1, col.names = c("resid", "rmsf"))
  x[, system := s]
  x
}
rmsf_all <- rbindlist(lapply(sys_names, load_rmsf))
rmsf_all <- merge(rmsf_all, TYPE_MAP, by = "system", all.x = TRUE)
rmsf_all[system == "APO", type := "APO"]
rmsf_all[system == "AZD5305", type := "Unk."]
rmsf_all[, domain := fifelse(resid %in% HD, "HD", fifelse(resid %in% ART, "ART", "other"))]
rmsf_all[, type := factor(type, levels = c("APO", "II", "III", "Unk."))]
rmsf_all[, system := factor(system, levels = sys_names)]

# (A) HD profiles
hd_prof <- rmsf_all[domain == "HD"]
pA <- ggplot(hd_prof, aes(resid, rmsf, color = system)) +
  geom_line(linewidth = 0.45) +
  scale_color_manual(values = cols) +
  labs(x = "HD residue (PARP1)", y = "RMSF (\u00C5)", title = "A  HD per-residue RMSF") +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9))

# (B) domain means
dom_mean <- rmsf_all[domain %in% c("HD", "ART"),
                     .(mean = mean(rmsf), sd = sd(rmsf)), by = .(system, type, domain)]
dom_mean[, domain := factor(domain, levels = c("HD", "ART"))]
pB <- ggplot(dom_mean, aes(system, mean, fill = domain)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.8), width = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = c(HD = "#C23531", ART = "#3D6BA8")) +
  labs(x = NULL, y = "Mean RMSF (\u00C5)", fill = "Domain",
       title = "B  Subdomain mean RMSF") +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9),
        axis.text.x = element_text(angle = 30, hjust = 1))

# (C) Type II - Type III dRMSF (HD domain)
dr <- rmsf_all[domain == "HD" & type %in% c("II", "III"),
               .(m = mean(rmsf)), by = .(resid, type)]
dr <- dcast(dr, resid ~ type, value.var = "m")
dr[, drmsf := II - III]
pC <- ggplot(dr, aes(resid, drmsf)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, color = "grey50") +
  geom_ribbon(aes(ymin = pmin(drmsf, 0), ymax = pmax(drmsf, 0)),
              fill = "#C23531", alpha = 0.35) +
  geom_line(linewidth = 0.5, color = "#9D2933") +
  labs(x = "HD residue (PARP1)", y = "\u0394RMSF (II \u2212 III) (\u00C5)",
       title = "C  Type II vs III \u0394RMSF") +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 9))

n_pos <- dr[drmsf > 0.3, .N]
cat(sprintf("Type II more flexible (dRMSF>0.3): %d HD residues\n", n_pos))
cat(sprintf("Type III more rigid side (dRMSF<-0.3): %d HD residues\n", dr[drmsf < -0.3, .N]))
top5 <- dr[order(-drmsf)][1:5]
cat("Top differentiated residues:\n"); print(top5)

p <- pA / (pB | pC) + plot_layout(heights = c(1, 1))
cairo_pdf("results/figures/Fig_S2_RMSF_All.pdf", width = 7.2, height = 6.0, pointsize = 8)
print(p)
dev.off()
cat("Saved: results/figures/Fig_S2_RMSF_All.pdf\n")
