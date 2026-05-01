# 03_diablo_model.R
# Tune and fit the DIABLO model. Tuning works out how many features
# to select per block per component — we try a range and let cross-
# validation pick the best. Then we fit the final model and check
# how well it actually classifies the three groups.
#
# Expect this to take 5-15 minutes depending on your machine.

library(mixOmics)
library(BiocParallel)

load("results/02_preprocessed_matrices.RData")

# ── design matrix ─────────────────────────────────────────────────────
# 0.1 = weak correlation prior between blocks — good for discovery
# we want the blocks to talk to each other a little but not be forced
design <- matrix(0.1,
                 nrow = 2, ncol = 2,
                 dimnames = list(c("proteomics", "metabolomics"),
                                 c("proteomics", "metabolomics")))
diag(design) <- 0

X <- list(proteomics = X_prot, metabolomics = X_met)

# ── tune keepX — how many features to select per component ───────────
set.seed(42)
cat("Tuning DIABLO — this takes a few minutes...\n")

tune_diablo <- tune.block.splsda(
  X          = X,
  Y          = Y,
  ncomp      = 2,
  test.keepX = list(
    proteomics   = c(5, 10, 15, 20, 25),
    metabolomics = c(5, 10, 15, 20, 25)
  ),
  design      = design,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 10,
  progressBar = TRUE
)

cat("Optimal keepX:\n")
print(tune_diablo$choice.keepX)

# ── fit the final model with tuned parameters ─────────────────────────
keepX_opt <- tune_diablo$choice.keepX

diablo_model <- block.splsda(
  X      = X,
  Y      = Y,
  ncomp  = 2,
  keepX  = keepX_opt,
  design = design
)

cat("Model fitted —",
    keepX_opt$proteomics[1], "proteins and",
    keepX_opt$metabolomics[1], "metabolites on comp 1\n")

# ── cross-validated performance check ────────────────────────────────
set.seed(42)
cat("Running performance evaluation...\n")

perf_diablo <- perf(
  diablo_model,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 10,
  progressBar = TRUE
)

cat("Error rates (Majority Vote):\n")
print(perf_diablo$MajorityVote.error.rate)

# ── pull out the selected biomarkers ─────────────────────────────────
prot_sel_c1 <- selectVar(diablo_model, block = "proteomics",   comp = 1)$proteomics$name
met_sel_c1  <- selectVar(diablo_model, block = "metabolomics", comp = 1)$metabolomics$name

cat("\nSelected proteins (comp 1):\n");   print(prot_sel_c1)
cat("\nSelected metabolites (comp 1):\n"); print(met_sel_c1)

# ── save everything ───────────────────────────────────────────────────
save(diablo_model, perf_diablo, tune_diablo, X, Y, matched_ids,
     keepX_opt, prot_sel_c1, met_sel_c1,
     file = "results/03_diablo_model.RData")
cat("Saved: results/03_diablo_model.RData\n")
