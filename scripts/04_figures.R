# 04_figures.R
# Generate the two main figures for the review article.
#
# Figure 1 — three-panel: sample plots for each block + arrow plot
#            showing how well the two blocks agree on sample positions
#
# Figure 2 — circos plot showing cross-omic correlations between
#            selected proteins and metabolites at r >= 0.6

library(mixOmics)

load("results/03_diablo_model.RData")

# colours: NCI = blue, MCI = orange, AD = red
group_cols <- c("#4575B4", "#74ADD1", "#D73027")

# ════════════════════════════════════════════════════════════════════
#  Figure 1 — sample plots + arrow plot
# ════════════════════════════════════════════════════════════════════

png("figures/Figure1_sample_plots.png",
    width = 3200, height = 1600, res = 300)

par(mfrow = c(1, 3))

# panel A — where do samples sit in proteomics space
plotIndiv(
  diablo_model,
  blocks        = "proteomics",
  comp          = c(1, 2),
  group         = Y,
  ind.names     = FALSE,
  legend        = TRUE,
  title         = "A  Proteomics",
  col.per.group = group_cols,
  pch           = 16,
  cex           = 1.2,
  ellipse       = TRUE
)

# panel B — same thing but in metabolomics space
plotIndiv(
  diablo_model,
  blocks        = "metabolomics",
  comp          = c(1, 2),
  group         = Y,
  ind.names     = FALSE,
  legend        = TRUE,
  title         = "B  Metabolomics",
  col.per.group = group_cols,
  pch           = 16,
  cex           = 1.2,
  ellipse       = TRUE
)

# panel C — arrows go from proteomics position to metabolomics position
# short arrows = the two blocks agree on where a sample sits
plotArrow(
  diablo_model,
  group         = Y,
  ind.names     = FALSE,
  legend        = TRUE,
  title         = "C  Block agreement",
  col.per.group = group_cols
)

dev.off()
cat("Figure 1 saved: figures/Figure1_sample_plots.png\n")

# ════════════════════════════════════════════════════════════════════
#  Figure 2 — circos plot
# ════════════════════════════════════════════════════════════════════
# cutoff 0.6 keeps only the strongest cross-omic correlations
# red lines = positive correlation, blue = negative

png("figures/Figure2_circos.png",
    width = 4500, height = 4500, res = 450)

par(mar = c(2, 2, 2, 2))

circosPlot(
  diablo_model,
  cutoff         = 0.6,
  line           = TRUE,
  color.blocks   = c("#4575B4", "#D73027"),
  color.cor      = c("#B2182B", "#2166AC"),
  size.labels    = 0.9,
  size.variables = 0.45,
  size.cex       = 0.8,
  var.adj        = 0.35
)

dev.off()
cat("Figure 2 saved: figures/Figure2_circos.png\n")

# ── also save the biomarker table ─────────────────────────────────────
# useful reference for the results section
prot_vars <- selectVar(diablo_model, block = "proteomics",   comp = 1)$proteomics
met_vars  <- selectVar(diablo_model, block = "metabolomics", comp = 1)$metabolomics

biomarker_table <- rbind(
  data.frame(feature   = rownames(prot_vars$value),
             loading   = round(prot_vars$value$value.var, 4),
             block     = "proteomics",
             component = 1),
  data.frame(feature   = rownames(met_vars$value),
             loading   = round(met_vars$value$value.var, 4),
             block     = "metabolomics",
             component = 1)
)

write.csv(biomarker_table, "results/biomarker_summary.csv", row.names = FALSE)
cat("Biomarker table saved: results/biomarker_summary.csv\n")
