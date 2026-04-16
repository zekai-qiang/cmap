# =============================================================================
# 03_longitudinal_lmm.R
#
# Linear mixed-effects model (LMM) of CMAP amplitude decline over time in
# ALS patients across four muscles: ADM, APB, BB, TA.
#
# Data are first normalised to baseline (percentage change from month 0).
# A random-slopes, fixed-intercept LMM is then fitted:
#
#   pct_change ~ time_months + (0 + time_months | study)
#
# This model allows each study to have its own rate of decline (random slope)
# while anchoring all studies at 0% change at baseline (fixed intercept).
# The model is fitted iteratively up to each timepoint (3, 6, 9, 12 months)
# to derive pooled estimates of decline at each observation window.
#
# An overlapping ggplot shows study-level trajectories alongside the pooled
# LMM estimate for each muscle, with Chan 2022 highlighted in red.
#
# Corresponds to: McKinnon, Qiang et al. (2024), Figure 2.
#
# Input:  data/longitudinal/lmm_<muscle>.csv
#         Expected columns: Study, X0, X3, X6, X9, X12
#         (one row per study; X0–X12 are mean CMAP amplitudes in mV)
#
# Output: LMM summary data frames (APB/ADM/BB/TA) and ggplot figures.
# =============================================================================

# --- Dependencies ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(ggeffects)
library(emmeans)
library(AICcmodavg)

# --- Data paths --------------------------------------------------------------
data_dir <- file.path("data", "longitudinal")

# --- Load per-muscle LMM data ------------------------------------------------
lmm_apb <- read.csv(file.path(data_dir, "lmm_apb.csv"))
lmm_adm <- read.csv(file.path(data_dir, "lmm_adm.csv"))
lmm_bb  <- read.csv(file.path(data_dir, "lmm_bb.csv"))
lmm_ta  <- read.csv(file.path(data_dir, "lmm_ta.csv"))

# --- Helper: convert to % change from baseline (month 0) --------------------
compute_pct_change <- function(df) {
  pct_df <- data.frame(
    sapply(
      df[, 2:ncol(df)],
      function(col) 100 * ((col - df$X0) / df$X0)
    )
  )
  pct_df$Study <- df$Study
  pct_df
}

# --- Helper: pivot to long format for LMM fitting ----------------------------
pivot_to_long <- function(pct_df) {
  long_df <- pct_df %>%
    pivot_longer(
      cols            = X0:X12,
      names_to        = "time_months",
      values_to       = "pct_change",
      values_transform = as.numeric
    )
  # Replace string column names (e.g. "X3") with integer month values
  n_timepoints  <- n_distinct(long_df$time_months)
  long_df$time_months <- rep(0:12, nrow(long_df) / 13)
  long_df
}

# --- Helper: fit LMM up to a given endpoint and return pooled estimate -------
# t:         endpoint in months (e.g. 3, 6, 9, 12)
# long_df:   long-format data frame for one muscle
# Returns:   c(pooled_mean, pooled_se, pooled_sd) at timepoint t
fit_lmm_at_endpoint <- function(t, long_df) {
  model <- lmer(
    pct_change ~ time_months + (0 + time_months | Study),
    data = filter(long_df, between(time_months, 0, t))
  )
  study_predictions <- crossing(
    Study       = factor(unique(long_df$Study)),
    time_months = 0:t
  ) %>%
    mutate(pct_change = predict(model, newdata = .))

  endpoint_preds <- filter(study_predictions, time_months == t)$pct_change
  c(
    pooled_mean = mean(endpoint_preds),
    pooled_se   = sd(endpoint_preds) / sqrt(length(endpoint_preds)),
    pooled_sd   = sd(endpoint_preds)
  )
}

# --- Helper: build LMM summary data frame for one muscle --------------------
build_lmm_summary <- function(long_df) {
  endpoint_months <- seq(3, 12, by = 3)
  estimates       <- sapply(endpoint_months, fit_lmm_at_endpoint, long_df = long_df)
  data.frame(
    time_months  = seq(0, 12, by = 3),
    pct_change   = c(0, estimates["pooled_mean", ]),
    std_error    = c(0, estimates["pooled_se",   ]),
    std_dev      = c(0, estimates["pooled_sd",   ])
  )
}

# --- Helper: plot LMM estimate with individual study trajectories ------------
# lmm_summary:  data frame from build_lmm_summary()
# long_df:      long-format study data for the same muscle
# muscle_label: character string for plot title (e.g. "APB")
plot_lmm_with_studies <- function(lmm_summary, long_df, muscle_label) {
  ggplot() +
    geom_point(
      data  = long_df,
      aes(x = time_months, y = pct_change, colour = Study)
    ) +
    geom_line(
      data  = filter(long_df, Study == "Chan 2022") %>% drop_na(),
      aes(x = time_months, y = pct_change, group = Study),
      colour    = "#c0392b",
      linewidth = 0.8
    ) +
    geom_point(
      data   = filter(long_df, Study == "Chan 2022") %>% drop_na(),
      aes(x  = time_months, y = pct_change),
      colour = "#c0392b",
      size   = 2.5
    ) +
    geom_line(
      data      = lmm_summary,
      aes(x     = time_months, y = pct_change),
      linewidth = 1,
      colour    = "black"
    ) +
    geom_point(
      data   = lmm_summary,
      aes(x  = time_months, y = pct_change),
      size   = 3,
      colour = "black"
    ) +
    labs(
      x     = "Time (months)",
      y     = "CMAP change from baseline (%)",
      title = sprintf("LMM-estimated CMAP decline — %s", muscle_label),
      caption = "Black: pooled LMM estimate. Red: Chan 2022. Colour: individual studies."
    ) +
    scale_x_continuous(breaks = seq(0, 12, by = 3)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
}

# =============================================================================
# Run LMM pipeline for all four muscles
# =============================================================================

# --- APB ---------------------------------------------------------------------
pct_change_apb  <- compute_pct_change(lmm_apb)
long_apb        <- pivot_to_long(pct_change_apb)
lmm_summary_apb <- build_lmm_summary(long_apb)
print(lmm_summary_apb)
plot_lmm_with_studies(lmm_summary_apb, long_apb, "APB")

# --- ADM ---------------------------------------------------------------------
pct_change_adm  <- compute_pct_change(lmm_adm)
long_adm        <- pivot_to_long(pct_change_adm)
lmm_summary_adm <- build_lmm_summary(long_adm)
print(lmm_summary_adm)
plot_lmm_with_studies(lmm_summary_adm, long_adm, "ADM")

# --- BB ----------------------------------------------------------------------
pct_change_bb  <- compute_pct_change(lmm_bb)
long_bb        <- pivot_to_long(pct_change_bb)
lmm_summary_bb <- build_lmm_summary(long_bb)
print(lmm_summary_bb)
plot_lmm_with_studies(lmm_summary_bb, long_bb, "BB")

# --- TA ----------------------------------------------------------------------
pct_change_ta  <- compute_pct_change(lmm_ta)
long_ta        <- pivot_to_long(pct_change_ta)
lmm_summary_ta <- build_lmm_summary(long_ta)
print(lmm_summary_ta)
plot_lmm_with_studies(lmm_summary_ta, long_ta, "TA")
