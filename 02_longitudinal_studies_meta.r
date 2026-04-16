# =============================================================================
# 02_longitudinal_studies_meta.R
#
# Longitudinal meta-analysis of CMAP amplitude decline in ALS patients across
# four muscles: ADM, APB, BB, TA. Pooled mean CMAP amplitudes are estimated at
# baseline and at 3, 6, 9, and 12 months using random-effects meta-analysis
# (REML, Knapp-Hartung adjustment).
#
# Studies reporting "NR" (not reported) at a given timepoint are excluded from
# that timepoint's analysis.
#
# Corresponds to: McKinnon, Qiang et al. (2024), Figure 2 and supplementary
# longitudinal figures.
#
# Input:  data/longitudinal/long_<muscle>.csv
#         Expected columns: study, n_base, baseline_mean, baseline_sd,
#                           n_3mo,  mean_3mo,  sd_3mo,
#                           n_6mo,  mean_6mo,  sd_6mo,
#                           n_9mo,  mean_9mo,  sd_9mo,
#                           n_12mo, mean_12mo, sd_12mo
#
# Output: Summary statistics and forest plots for each muscle × timepoint.
# =============================================================================

# --- Dependencies ------------------------------------------------------------
library(meta)
library(tidyverse)
library(metafor)

# --- Data paths --------------------------------------------------------------
data_dir <- file.path("data", "longitudinal")

# --- Load longitudinal data --------------------------------------------------
long_adm <- read_csv(file.path(data_dir, "long_adm.csv"), show_col_types = FALSE)
long_apb <- read_csv(file.path(data_dir, "long_apb.csv"), show_col_types = FALSE)
long_bb  <- read_csv(file.path(data_dir, "long_bb.csv"),  show_col_types = FALSE)
long_ta  <- read_csv(file.path(data_dir, "long_ta.csv"),  show_col_types = FALSE)

# --- Helper: random-effects meta-analysis for a single timepoint -------------
# df:       longitudinal data frame for one muscle
# n_col:    column name for sample size
# mean_col: column name for mean CMAP
# sd_col:   column name for SD
# Rows where n_col, mean_col, or sd_col == "NR" are automatically excluded.
run_longitudinal_meta <- function(df, n_col, mean_col, sd_col) {
  df %>%
    filter(
      .data[[mean_col]] != "NR",
      .data[[sd_col]]   != "NR",
      .data[[n_col]]    != "NR"
    ) %>%
    metamean(
      n          = as.numeric(.data[[n_col]]),
      mean       = as.numeric(.data[[mean_col]]),
      sd         = as.numeric(.data[[sd_col]]),
      sm         = "MRAW",
      studlab    = study,
      fixed      = FALSE,
      random     = TRUE,
      method.tau = "REML",
      hakn       = TRUE
    )
}

# --- Helper: run and report all timepoints for one muscle --------------------
analyse_muscle_longitudinal <- function(df, muscle_label) {
  timepoints <- list(
    baseline = list(n = "n_base",  mean = "baseline_mean", sd = "baseline_sd"),
    mo3      = list(n = "n_3mo",   mean = "mean_3mo",      sd = "sd_3mo"),
    mo6      = list(n = "n_6mo",   mean = "mean_6mo",      sd = "sd_6mo"),
    mo9      = list(n = "n_9mo",   mean = "mean_9mo",      sd = "sd_9mo"),
    mo12     = list(n = "n_12mo",  mean = "mean_12mo",     sd = "sd_12mo")
  )

  results <- list()
  for (tp_name in names(timepoints)) {
    tp     <- timepoints[[tp_name]]
    meta_result <- run_longitudinal_meta(df,
                                         n_col    = tp$n,
                                         mean_col = tp$mean,
                                         sd_col   = tp$sd)
    message(sprintf("\n--- %s | %s ---", muscle_label, tp_name))
    print(summary(meta_result))
    forest.meta(meta_result,
                main = sprintf("%s CMAP — %s", muscle_label, tp_name))
    results[[tp_name]] <- meta_result
  }
  invisible(results)
}

# =============================================================================
# Run longitudinal meta-analyses for all four muscles
# =============================================================================
adm_meta_results <- analyse_muscle_longitudinal(long_adm, "ADM")
apb_meta_results <- analyse_muscle_longitudinal(long_apb, "APB")
bb_meta_results  <- analyse_muscle_longitudinal(long_bb,  "BB")
ta_meta_results  <- analyse_muscle_longitudinal(long_ta,  "TA")
