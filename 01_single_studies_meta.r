# =============================================================================
# 01_single_studies_meta.R
#
# Cross-sectional meta-analysis of CMAP amplitudes in ALS patients and
# healthy controls (HCs) across four muscles: ADM, APB, FDI, TA.
#
# Each muscle is analysed separately for HC and ALS cohorts using a
# random-effects meta-analysis (restricted maximum likelihood, REML),
# with the Knapp-Hartung adjustment (hakn = TRUE). Forest plots are
# generated for each analysis.
#
# Corresponds to: McKinnon, Qiang et al. (2024), Figure 1 and related
# supplementary figures (human cross-sectional data).
#
# Input:  data/single_studies/<muscle>_hc.csv
#         data/single_studies/<muscle>_als.csv
# Output: Forest plots (printed to active graphics device)
# =============================================================================

# --- Dependencies ------------------------------------------------------------
library(meta)
library(tidyverse)
library(metafor)
library(grid)

# --- Data paths --------------------------------------------------------------
# Set this to the directory containing the data/ folder
data_dir <- file.path("data", "single_studies")

# --- Helper: load cross-sectional data pair ----------------------------------
load_cmap_data <- function(muscle_code) {
  hc  <- read_csv(file.path(data_dir, paste0(muscle_code, "_hc.csv")),
                  show_col_types = FALSE)
  als <- read_csv(file.path(data_dir, paste0(muscle_code, "_als.csv")),
                  show_col_types = FALSE)
  list(hc = hc, als = als)
}

# --- Helper: run random-effects meta-analysis --------------------------------
# Expects data frame with columns: study, n, cmap_mean, cmap_sd
run_meta_mean <- function(df, n_col, mean_col, sd_col) {
  df %>%
    metamean(
      n          = .data[[n_col]],
      mean       = .data[[mean_col]],
      sd         = .data[[sd_col]],
      sm         = "MRAW",
      studlab    = study,
      fixed      = FALSE,
      random     = TRUE,
      method.tau = "REML",
      hakn       = TRUE
    )
}

# --- Helper: forest plot with title ------------------------------------------
plot_forest <- function(meta_obj, title_text) {
  forest.meta(meta_obj, sortvar = TE)
  grid.text(title_text, x = 0.5, y = 0.93, gp = gpar(cex = 1))
}

# =============================================================================
# 1. ADM — Abductor Digiti Minimi
# =============================================================================
adm_data <- load_cmap_data("adm")

adm_hc_meta <- run_meta_mean(adm_data$hc,
                              n_col    = "n_hc",
                              mean_col = "cmap_hc_mean",
                              sd_col   = "cmap_hc_sd")
summary(adm_hc_meta)
plot_forest(adm_hc_meta, "Pooled mean CMAP amplitude of ADM in healthy controls (mV)")

adm_als_meta <- run_meta_mean(adm_data$als,
                               n_col    = "n_als",
                               mean_col = "cmap_als_mean",
                               sd_col   = "cmap_als_sd")
summary(adm_als_meta)
plot_forest(adm_als_meta, "Pooled mean CMAP amplitude of ADM in ALS (mV)")

# =============================================================================
# 2. APB — Abductor Pollicis Brevis
# =============================================================================
apb_data <- load_cmap_data("apb")

apb_hc_meta <- run_meta_mean(apb_data$hc,
                              n_col    = "n_hc",
                              mean_col = "cmap_hc_mean",
                              sd_col   = "cmap_hc_sd")
summary(apb_hc_meta)
plot_forest(apb_hc_meta, "Pooled mean CMAP amplitude of APB in healthy controls (mV)")

apb_als_meta <- run_meta_mean(apb_data$als,
                               n_col    = "n_als",
                               mean_col = "cmap_als_mean",
                               sd_col   = "cmap_als_sd")
summary(apb_als_meta)
plot_forest(apb_als_meta, "Pooled mean CMAP amplitude of APB in ALS (mV)")

# =============================================================================
# 3. FDI — First Dorsal Interosseous
# =============================================================================
fdi_data <- load_cmap_data("fdi")

fdi_hc_meta <- run_meta_mean(fdi_data$hc,
                              n_col    = "n_hc",
                              mean_col = "cmap_hc_mean",
                              sd_col   = "cmap_hc_sd")
summary(fdi_hc_meta)
plot_forest(fdi_hc_meta, "Pooled mean CMAP amplitude of FDI in healthy controls (mV)")

fdi_als_meta <- run_meta_mean(fdi_data$als,
                               n_col    = "n_als",
                               mean_col = "cmap_als_mean",
                               sd_col   = "cmap_als_sd")
summary(fdi_als_meta)
plot_forest(fdi_als_meta, "Pooled mean CMAP amplitude of FDI in ALS (mV)")

# =============================================================================
# 4. TA — Tibialis Anterior
# =============================================================================
ta_data <- load_cmap_data("ta")

ta_hc_meta <- run_meta_mean(ta_data$hc,
                             n_col    = "n_hc",
                             mean_col = "cmap_hc_mean",
                             sd_col   = "cmap_hc_sd")
summary(ta_hc_meta)
plot_forest(ta_hc_meta, "Pooled mean CMAP amplitude of TA in healthy controls (mV)")

ta_als_meta <- run_meta_mean(ta_data$als,
                              n_col    = "n_als",
                              mean_col = "cmap_als_mean",
                              sd_col   = "cmap_als_sd")
summary(ta_als_meta)
plot_forest(ta_als_meta, "Pooled mean CMAP amplitude of TA in ALS (mV)")
