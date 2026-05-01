# 02_preprocess.R
# Turn the raw data files into clean matrices that mixOmics can actually
# use. Proteomics is already log2 TMT so we leave it alone. Metabolomics
# gets a log2 transform because the raw concentrations are right-skewed.
# Missing values are filled with column medians — simple and good enough
# for this kind of exploratory integration.

library(readr)
library(dplyr)
library(mixOmics)

# dplyr::select gets masked by igraph when mixOmics loads, pin it now
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange

# ── paths ─────────────────────────────────────────────────────────────
base_prot <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Proteomics/"
base_met  <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Metabolomics/Metabolon/Data/"

# load matched IDs from step 1
load("results/01_matched_samples.RData")

prot <- read_csv(paste0(base_prot, "ROSMAP_proteomics_merged_clean.csv"), show_col_types = FALSE)
met  <- read_csv(paste0(base_met,  "ROSMAP_Metabolon_HD4_Brain514_assay_data.csv"),
                 show_col_types = FALSE, locale = locale(encoding = "latin1"))

# ── proteomics matrix ─────────────────────────────────────────────────
# drop the ID/metadata columns, keep only the numeric protein columns
prot_meta_cols <- c("SampleID", "batch.channel", "projid", "IndividualID",
                    "SpecimenID", "MulticonsensusStudyFileID")

X_prot <- prot |>
  dplyr::filter(IndividualID %in% matched_ids) |>
  dplyr::arrange(match(IndividualID, matched_ids)) |>
  dplyr::select(-dplyr::all_of(intersect(prot_meta_cols, colnames(prot)))) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
  as.matrix()

rownames(X_prot) <- matched_ids

# DxNumeric snuck in as a column — it's a clinical label, not a protein
X_prot <- X_prot[, !colnames(X_prot) %in% "DxNumeric.Emory2017"]

# fill missing with column median — proteins are already log2 so no transform
col_medians_p <- apply(X_prot, 2, median, na.rm = TRUE)
for (j in seq_len(ncol(X_prot))) {
  X_prot[is.na(X_prot[, j]), j] <- col_medians_p[j]
}

# catch any Inf values that sneak through
X_prot[!is.finite(X_prot)] <- NA
col_medians_p2 <- apply(X_prot, 2, median, na.rm = TRUE)
for (j in seq_len(ncol(X_prot))) {
  X_prot[is.na(X_prot[, j]), j] <- col_medians_p2[j]
}

# drop near-zero variance proteins — mixOmics will refuse them anyway
nzv_prot <- nearZeroVar(X_prot)
if (length(nzv_prot$Position) > 0) X_prot <- X_prot[, -nzv_prot$Position]
colnames(X_prot) <- make.names(colnames(X_prot), unique = TRUE)

cat("Proteomics matrix:", nrow(X_prot), "x", ncol(X_prot), "\n")

# ── metabolomics matrix ───────────────────────────────────────────────
# character columns are metadata, everything else is a measurement
met_char_cols <- names(met)[sapply(met, is.character)]

X_met <- met |>
  dplyr::filter(individualID %in% matched_ids) |>
  dplyr::mutate(individualID = as.character(individualID)) |>
  dplyr::arrange(match(individualID, matched_ids)) |>
  dplyr::select(-dplyr::all_of(met_char_cols)) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
  as.matrix()

rownames(X_met) <- matched_ids

# log2 transform — zeros and negatives go to NA first
X_met[X_met <= 0] <- NA
X_met <- log2(X_met)

# impute missing with column median
col_medians_m <- apply(X_met, 2, median, na.rm = TRUE)
for (j in seq_len(ncol(X_met))) {
  X_met[is.na(X_met[, j]), j] <- col_medians_m[j]
}

X_met[!is.finite(X_met)] <- NA
col_medians_m2 <- apply(X_met, 2, median, na.rm = TRUE)
for (j in seq_len(ncol(X_met))) {
  X_met[is.na(X_met[, j]), j] <- col_medians_m2[j]
}

# same near-zero variance filter
nzv_met <- nearZeroVar(X_met)
if (length(nzv_met$Position) > 0) X_met <- X_met[, -nzv_met$Position]
colnames(X_met) <- make.names(colnames(X_met), unique = TRUE)

cat("Metabolomics matrix:", nrow(X_met), "x", ncol(X_met), "\n")

# ── decode metabolite IDs to readable names ───────────────────────────
met_dict <- read_csv(
  "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Metabolomics/Metabolon/QC/ROSMAP Metabolon HD4 Data Dictionary.csv",
  show_col_types = FALSE
)

# CHEM_ID matches the number after "X" in column names e.g. X100001620
met_id_map <- met_dict |>
  dplyr::mutate(col_name = paste0("X", CHEM_ID)) |>
  dplyr::select(col_name, CHEMICAL_NAME)

col_lookup    <- setNames(met_id_map$CHEMICAL_NAME, met_id_map$col_name)
new_met_names <- ifelse(
  colnames(X_met) %in% names(col_lookup) & !is.na(col_lookup[colnames(X_met)]),
  col_lookup[colnames(X_met)],
  colnames(X_met)
)
colnames(X_met) <- make.unique(new_met_names, sep = "_v")

# quick sanity check before handing off to DIABLO
cat("Row order matches:", all(rownames(X_prot) == rownames(X_met)), "\n")
cat("Missing in prot:", sum(is.na(X_prot)), "| Missing in met:", sum(is.na(X_met)), "\n")

save(X_prot, X_met, Y, matched_ids,
     file = "results/02_preprocessed_matrices.RData")
cat("Saved: results/02_preprocessed_matrices.RData\n")
