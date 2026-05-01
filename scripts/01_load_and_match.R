# 01_load_and_match.R
# Pull in the proteomics, metabolomics and clinical files, then figure out
# which samples we actually have across all three. That overlap is our
# working dataset for the rest of the analysis.
#
# Blocks used : proteomics + metabolomics
# Outcome     : cognitive diagnosis (NCI / MCI / AD) from cogdx column
# Note        : GEX and variant blocks were excluded — their sample ID
#               mapping files were not available in this data export.

library(readr)
library(dplyr)

# ── paths — update these if you move the data folder ─────────────────
base_var  <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Gene Variants Dataset and Analysis/"
base_prot <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Proteomics/"
base_met  <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Metabolomics/Metabolon/Data/"

# ── load the three files we need ─────────────────────────────────────
clin <- read_csv(paste0(base_var,  "ROSMAP_clinical.csv"),              show_col_types = FALSE)
prot <- read_csv(paste0(base_prot, "ROSMAP_proteomics_merged_clean.csv"), show_col_types = FALSE)
met  <- read_csv(paste0(base_met,  "ROSMAP_Metabolon_HD4_Brain514_assay_data.csv"),
                 show_col_types = FALSE, locale = locale(encoding = "latin1"))

# ── build the outcome table — only keeping the three clean cogdx groups
# cogdx 1 = no impairment, 2 = MCI, 4 = AD; 3/5/6 are too small to use
id_bridge <- clin |>
  dplyr::select(projid, individualID, cogdx) |>
  dplyr::mutate(
    projid       = as.character(projid),
    individualID = as.character(individualID),
    group = dplyr::case_when(
      cogdx == 1 ~ "NCI",
      cogdx == 2 ~ "MCI",
      cogdx == 4 ~ "AD",
      TRUE       ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(individualID), !is.na(group))

# ── pull the sample IDs that exist in each block ─────────────────────
prot_ids <- prot |>
  dplyr::filter(!is.na(IndividualID)) |>
  dplyr::distinct(IndividualID) |>
  dplyr::pull(IndividualID)

met_ids <- met |>
  dplyr::filter(!is.na(individualID)) |>
  dplyr::mutate(individualID = as.character(individualID)) |>
  dplyr::distinct(individualID) |>
  dplyr::pull(individualID)

# keep only samples that appear in all three — prot, met, and clinical
matched_ids <- Reduce(intersect, list(prot_ids, met_ids, id_bridge$individualID))
cat("Matched samples across all blocks:", length(matched_ids), "\n")

# outcome vector locked to the same order as matched_ids
outcome_df <- id_bridge |>
  dplyr::filter(individualID %in% matched_ids) |>
  dplyr::arrange(match(individualID, matched_ids))

Y <- factor(outcome_df$group, levels = c("NCI", "MCI", "AD"))
cat("Group breakdown:\n")
print(table(Y))

# save the matched IDs and outcome so the next script can pick them up
save(matched_ids, Y, outcome_df,
     file = "results/01_matched_samples.RData")
cat("Saved: results/01_matched_samples.RData\n")
