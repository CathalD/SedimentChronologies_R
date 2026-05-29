# =============================================================================
# RERCA | Script 05: Multi-Model Summary and Comparison
# =============================================================================
#
# OVERVIEW
# --------
# This script loads the saved outputs from scripts 01–04 and produces:
#
#   1. A combined age-depth plot overlaying all models
#   2. An activity profile plot (Pb-210 vs depth)
#   3. A sedimentation rate comparison across models
#   4. A summary statistics table (CSV)
#   5. An age comparison table at common depth intervals (CSV)
#   6. A console report with statistical context for each metric
#
# All outputs are saved to RERCA/output/
#
# PREREQUISITES
# -------------
# Run scripts 01–04 first so their output CSVs exist in RERCA/output/.
# Script 04 (Bayesian/Stan) is optional — the summary will work without it.
#
# STEPS
# -----
# 1. Install / load packages
# 2. Load all model outputs
# 3. Combined age-depth plot
# 4. Activity profile
# 5. Sedimentation rate profile
# 6. Age comparison table at common depths
# 7. Model statistics table
# 8. Console statistical report
# =============================================================================


# =============================================================================
# STEP 1 — Install and load packages
# =============================================================================

if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("scales",  quietly = TRUE)) install.packages("scales")

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

output_dir <- "RERCA/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# STEP 2 — Load all model outputs
# =============================================================================

load_if_exists <- function(path) {
  if (file.exists(path)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    cat("Loaded:", path, sprintf("(%d rows)\n", nrow(df)))
    df
  } else {
    cat("Not found (skip):", path, "\n")
    NULL
  }
}

cat("Loading model outputs from", output_dir, "\n\n")

crs_ages     <- load_if_exists(file.path(output_dir, "01_crs_ages.csv"))
rplum_ages   <- load_if_exists(file.path(output_dir, "02_rplum_ages.csv"))
serac_ages   <- load_if_exists(file.path(output_dir, "03_serac_ages.csv"))
bayes_ages   <- load_if_exists(file.path(output_dir, "04_bayesian_ages.csv"))
activity_df  <- load_if_exists(file.path(output_dir, "01_activity_profile.csv"))

crs_stats    <- load_if_exists(file.path(output_dir, "01_crs_stats.csv"))
rplum_stats  <- load_if_exists(file.path(output_dir, "02_rplum_stats.csv"))
serac_stats  <- load_if_exists(file.path(output_dir, "03_serac_stats.csv"))
bayes_stats  <- load_if_exists(file.path(output_dir, "04_bayesian_stats.csv"))

# Harmonise column names across all age outputs and stack into one data frame
harmonise <- function(df, model_label) {
  if (is.null(df)) return(NULL)
  # Keep only the columns we need; rename if necessary
  depth_col <- intersect(c("depth_cm", "depth"), names(df))[1]
  df %>%
    rename_with(~ "depth_cm", any_of(c("depth", "depth_cm"))) %>%
    select(depth_cm,
           age        = any_of(c("age", "mean")),
           age_min_yr = any_of(c("age_min_yr", "min.95")),
           age_max_yr = any_of(c("age_max_yr", "max.95"))) %>%
    mutate(model = model_label)
}

all_ages <- bind_rows(
  harmonise(crs_ages,   "pb210 CRS"),
  harmonise(rplum_ages, "rplum"),
  serac_ages %>%
    rename(depth_cm = any_of(c("depth_cm","depth"))) %>%
    select(depth_cm, age, age_min_yr, age_max_yr, model) %>%
    mutate(model = paste0("serac ", sub("serac_", "", model))),
  harmonise(bayes_ages, "Bayesian (Stan)")
) %>%
  filter(!is.na(age), !is.na(depth_cm), age >= 0)

n_models <- length(unique(all_ages$model))
cat(sprintf("\n%d model(s) loaded: %s\n\n",
            n_models, paste(unique(all_ages$model), collapse = ", ")))

if (n_models == 0) stop("No model outputs found. Run scripts 01–04 first.")

# Colour palette — up to 6 models
model_colours <- c(
  "pb210 CRS"      = "#2c7bb6",
  "rplum"          = "#d7191c",
  "serac CRS"      = "#756bb1",
  "serac CFCS"     = "#1a9641",
  "Bayesian (Stan)"= "#f16913"
)
# Assign colours to whatever models are present
present_models  <- unique(all_ages$model)
palette_colours <- model_colours[names(model_colours) %in% present_models]
missing_models  <- setdiff(present_models, names(model_colours))
if (length(missing_models) > 0) {
  extra <- setNames(hue_pal()(length(missing_models)), missing_models)
  palette_colours <- c(palette_colours, extra)
}


# =============================================================================
# STEP 3 — Combined age-depth plot
# =============================================================================

p_all_ages <- ggplot(all_ages, aes(x = age, y = depth_cm, colour = model, fill = model)) +
  geom_ribbon(aes(xmin = age_min_yr, xmax = age_max_yr),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Age (years before coring)") +
  scale_colour_manual(values = palette_colours) +
  scale_fill_manual(values = palette_colours, guide = "none") +
  labs(title    = "Pb-210 Age-Depth Model Comparison",
       subtitle = "Shaded bands = 95% uncertainty intervals",
       colour   = "Model") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        legend.title     = element_text(face = "bold"))

print(p_all_ages)
ggsave(file.path(output_dir, "05_age_depth_comparison.pdf"),
       p_all_ages, width = 7, height = 6)
cat("Saved: 05_age_depth_comparison.pdf\n")


# =============================================================================
# STEP 4 — Pb-210 Activity Profile
# =============================================================================

if (!is.null(activity_df)) {
  bg_line <- unique(activity_df$pb210_supported)[1]

  act_long <- activity_df %>%
    select(depth_cm, total = pb210_total_dpm_g, unsupported = pb210_unsupported) %>%
    pivot_longer(cols = c(total, unsupported),
                 names_to = "fraction", values_to = "activity") %>%
    left_join(
      activity_df %>%
        select(depth_cm,
               total_err = pb210_error_dpm_g,
               unsup_err = pb210_unsupported_err),
      by = "depth_cm"
    ) %>%
    mutate(err = ifelse(fraction == "total", total_err, unsup_err),
           fraction = recode(fraction,
                             total       = "Total Pb-210",
                             unsupported = "Unsupported Pb-210"))

  p_activity <- ggplot(act_long, aes(x = activity, y = depth_cm, colour = fraction)) +
    geom_point(size = 2.5) +
    geom_errorbarh(aes(xmin = activity - err, xmax = activity + err), height = 0.3) +
    geom_vline(xintercept = bg_line, linetype = "dashed", colour = "grey40") +
    annotate("text", x = bg_line, y = max(activity_df$depth_cm) * 0.9,
             label = sprintf("Background\n%.3f DPM/g", bg_line),
             hjust = -0.1, size = 3, colour = "grey40") +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "Activity (DPM / g dry wt)") +
    scale_colour_manual(values = c("Total Pb-210"       = "#2c7bb6",
                                   "Unsupported Pb-210" = "#d7191c")) +
    labs(title    = "Pb-210 Activity Profile",
         subtitle = "Dashed line = background (supported) activity",
         colour   = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")

  print(p_activity)
  ggsave(file.path(output_dir, "05_activity_profile.pdf"),
         p_activity, width = 6, height = 5)
  cat("Saved: 05_activity_profile.pdf\n")
}


# =============================================================================
# STEP 5 — Sedimentation rate profiles
# =============================================================================
# Only models that report sedimentation rate at each section (CRS from script 01)
# have a per-depth rate; serac's CFCS is a single value.

if (!is.null(crs_ages) && "sed_rate_cm_yr" %in% names(crs_ages)) {
  sed_df <- crs_ages %>%
    filter(!is.na(sed_rate_cm_yr), sed_rate_cm_yr > 0) %>%
    mutate(model = "pb210 CRS")

  p_sed <- ggplot(sed_df, aes(x = sed_rate_cm_yr, y = depth_cm)) +
    geom_step(orientation = "y", colour = palette_colours["pb210 CRS"], linewidth = 1) +
    geom_point(colour = palette_colours["pb210 CRS"], size = 2.5) +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "Sedimentation rate (cm / yr)") +
    labs(title    = "CRS Sedimentation Rate Profile",
         subtitle = "Variable rate allowed under the CRS assumption") +
    theme_bw(base_size = 13)
  print(p_sed)
  ggsave(file.path(output_dir, "05_sedimentation_rate.pdf"),
         p_sed, width = 6, height = 5)
  cat("Saved: 05_sedimentation_rate.pdf\n")
}


# =============================================================================
# STEP 6 — Age comparison table at common depth intervals
# =============================================================================
# Interpolate each model's age to a regular 2-cm depth grid for direct comparison.

depth_grid <- seq(0, floor(max(all_ages$depth_cm, na.rm = TRUE)), by = 2)

interp_model <- function(df, model_name, depths) {
  df_m <- df %>% filter(model == model_name) %>% arrange(depth_cm)
  if (nrow(df_m) < 2) return(NULL)
  data.frame(
    depth_cm  = depths,
    age       = approx(df_m$depth_cm, df_m$age,        xout = depths, rule = 1)$y,
    age_lo    = approx(df_m$depth_cm, df_m$age_min_yr, xout = depths, rule = 1)$y,
    age_hi    = approx(df_m$depth_cm, df_m$age_max_yr, xout = depths, rule = 1)$y,
    model     = model_name
  )
}

age_grid <- bind_rows(lapply(unique(all_ages$model),
                             function(m) interp_model(all_ages, m, depth_grid)))

age_wide <- age_grid %>%
  select(depth_cm, model, age) %>%
  pivot_wider(names_from = model, values_from = age,
              names_glue = "{model}_age_yr") %>%
  mutate(across(where(is.numeric), \(x) round(x, 1)))

cat("\n=== Age Comparison at 2-cm Depth Intervals ===\n")
print(age_wide)
write.csv(age_wide, file.path(output_dir, "05_age_comparison_table.csv"), row.names = FALSE)
cat("Saved: 05_age_comparison_table.csv\n")

# Also compute inter-model age range at each depth (measure of agreement)
age_wide <- age_wide %>%
  rowwise() %>%
  mutate(
    model_age_min  = min(c_across(ends_with("_age_yr")), na.rm = TRUE),
    model_age_max  = max(c_across(ends_with("_age_yr")), na.rm = TRUE),
    model_age_range = model_age_max - model_age_min
  ) %>%
  ungroup()

cat("\nInter-model age spread at each depth (smaller = better agreement):\n")
print(age_wide %>%
  select(depth_cm, model_age_min, model_age_max, model_age_range) %>%
  mutate(across(where(is.numeric), \(x) round(x, 1))))


# =============================================================================
# STEP 7 — Model statistics table
# =============================================================================

stats_all <- bind_rows(crs_stats, rplum_stats, serac_stats, bayes_stats)

if (nrow(stats_all) > 0) {
  cat("\n=== Model Statistics Summary ===\n")
  print(stats_all)
  write.csv(stats_all, file.path(output_dir, "05_model_statistics.csv"), row.names = FALSE)
  cat("Saved: 05_model_statistics.csv\n")
}


# =============================================================================
# STEP 8 — Console statistical report
# =============================================================================

cat("\n")
cat("==========================================================================\n")
cat("  RERCA Pb-210 CHRONOLOGY — STATISTICAL REPORT\n")
cat("==========================================================================\n\n")

# --- Activity profile assessment ---
if (!is.null(activity_df)) {
  n_measured      <- nrow(activity_df)
  max_depth_dated <- max(all_ages$depth_cm, na.rm = TRUE)
  bg              <- unique(activity_df$pb210_supported)[1]
  surface_act     <- activity_df$pb210_total_dpm_g[which.min(activity_df$depth_cm)]
  dynamic_range   <- surface_act / bg

  cat("ACTIVITY PROFILE\n")
  cat("----------------\n")
  cat(sprintf("  Measured sections:          %d\n", n_measured))
  cat(sprintf("  Background (supported):     %.3f DPM/g\n", bg))
  cat(sprintf("  Surface total activity:     %.3f DPM/g\n", surface_act))
  cat(sprintf("  Dynamic range (surf/bg):    %.1fx\n", dynamic_range))
  cat(sprintf("  Datable depth horizon:      %.1f cm\n", max_depth_dated))
  cat("\n")
  cat("  What this means:\n")
  cat("  The dynamic range is the ratio of surface to background activity.\n")
  cat("  A higher ratio gives better model constraints.\n")
  cat("  Ratios < 3 indicate limited signal; > 5 is robust for CRS.\n\n")
}

# --- Datable horizon ---
cat("DATABLE HORIZON\n")
cat("---------------\n")
for (m in unique(all_ages$model)) {
  sub <- all_ages %>% filter(model == m)
  oldest <- sub %>% filter(age == max(age, na.rm = TRUE))
  cat(sprintf("  %-20s  %.1f cm  |  %.0f yr BP  (95%% CI: %.0f – %.0f yr)\n",
              m,
              max(sub$depth_cm, na.rm = TRUE),
              oldest$age[1],
              oldest$age_min_yr[1],
              oldest$age_max_yr[1]))
}
cat("\n  What this means:\n")
cat("  The datable horizon is the deepest depth where unsupported Pb-210\n")
cat("  is distinguishable from background. Sections below this depth\n")
cat("  cannot be reliably dated with Pb-210 and should not be used\n")
cat("  for RERCA calculations.\n\n")

# --- Inter-model agreement ---
cat("INTER-MODEL AGREEMENT\n")
cat("---------------------\n")
if (nrow(age_wide) > 0) {
  mean_spread <- mean(age_wide$model_age_range, na.rm = TRUE)
  max_spread  <- max(age_wide$model_age_range,  na.rm = TRUE)
  cat(sprintf("  Mean age spread across models:  %.1f yr\n", mean_spread))
  cat(sprintf("  Maximum age spread:             %.1f yr  (at %.0f cm)\n",
              max_spread,
              age_wide$depth_cm[which.max(age_wide$model_age_range)]))
  cat("\n  What this means:\n")
  cat("  Small spread (< 10 yr) = strong model agreement.\n")
  cat("  Large spread at a particular depth may indicate a\n")
  cat("  depositional change (e.g. flood layer, erosion) that\n")
  cat("  violates one or more model assumptions.\n\n")
}

# --- CFCS linearity ---
if (!is.null(serac_stats) && "cfcs_r2" %in% names(serac_stats)) {
  cfcs_row <- serac_stats %>% filter(!is.na(cfcs_r2))
  if (nrow(cfcs_row) > 0) {
    cat("CFCS MODEL (serac)\n")
    cat("------------------\n")
    cat(sprintf("  SAR:  %.3f mm/yr  (+/- %.3f)\n",
                abs(cfcs_row$cfcs_sar_mm_yr[1]),
                abs(cfcs_row$cfcs_sar_err[1])))
    cat(sprintf("  R²:   %.3f\n", cfcs_row$cfcs_r2[1]))
    cat("\n  What this means:\n")
    cat("  R² measures how well a single constant sedimentation rate fits\n")
    cat("  the log-linear Pb-210 decay profile.\n")
    cat("  R² > 0.9 suggests constant sedimentation — CFCS appropriate.\n")
    cat("  R² < 0.8 suggests variable rates — CRS model preferred.\n\n")
  }
}

# --- Bayesian diagnostics ---
if (!is.null(bayes_stats) && "n_divergent" %in% names(bayes_stats)) {
  cat("BAYESIAN MODEL (Stan)\n")
  cat("---------------------\n")
  cat(sprintf("  Divergent transitions:  %d\n", bayes_stats$n_divergent[1]))
  cat(sprintf("  phi (supported) mean:   %.3f DPM/g  (95%% CI: %.3f – %.3f)\n",
              bayes_stats$phi_mean_dpm_g[1],
              bayes_stats$phi_95_lo[1],
              bayes_stats$phi_95_hi[1]))
  cat("\n  What this means:\n")
  cat("  Divergent transitions > 10 indicate poor MCMC mixing — results\n")
  cat("  should be treated cautiously. Consider tightening priors.\n")
  cat("  0 divergences = sampler converged well.\n\n")
}

# --- RERCA usability verdict ---
cat("==========================================================================\n")
cat("  USABILITY FOR RERCA (Carbon Accumulation Rate)\n")
cat("==========================================================================\n\n")
cat("  To convert an age-depth model into a carbon accumulation rate (RERCA),\n")
cat("  you need:\n")
cat("    1. A reliable age-depth model (see above)\n")
cat("    2. Bulk density measurements for all dated sections\n")
cat("    3. LOI or % organic carbon for all dated sections\n\n")
cat("  Formula:\n")
cat("    RERCA (g C m-2 yr-1) = LOI_fraction * bulk_density (g/cm3)\n")
cat("                           * sedimentation_rate (cm/yr) * 10000\n\n")
cat("  Key requirements for a usable age-depth model:\n")
cat("    - At least 3 dated points above the datable horizon\n")
cat("    - 95% CI width < 50% of estimated age at each point\n")
cat("    - Consistent ages across at least 2 models\n")
cat("    - Surface age consistent with coring year\n\n")

# Auto-check key requirements
n_dated_pts <- all_ages %>% group_by(model) %>% summarise(n = n())
ci_width    <- all_ages %>% mutate(ci_pct = (age_max_yr - age_min_yr) / age * 100) %>%
               group_by(model) %>% summarise(mean_ci_pct = mean(ci_pct, na.rm = TRUE))

cat("  Auto-checks:\n")
for (m in unique(all_ages$model)) {
  nd  <- n_dated_pts$n[n_dated_pts$model == m]
  ci  <- ci_width$mean_ci_pct[ci_width$model == m]
  cat(sprintf("    %-20s  dated points: %2d  |  mean 95%% CI width: %.0f%% of age\n",
              m, nd, ci))
}
cat("\n")
cat("  Review the comparison plot and table, then select the model whose\n")
cat("  age-depth curve best matches your knowledge of the site history.\n\n")
cat("==========================================================================\n")

cat("\nAll outputs saved to:", output_dir, "\n")
cat("Script 05 complete.\n")


# =============================================================================
# STEP 9 — Build and render Quarto summary report
# =============================================================================
# All data objects are already in memory. We:
#   (a) save a snapshot to RData so the .qmd can load them cleanly,
#   (b) write the .qmd file,
#   (c) render it with quarto::quarto_render().

if (!requireNamespace("quarto", quietly = TRUE)) install.packages("quarto")

# --- 9a. Save report data snapshot ---
report_data_path <- file.path(output_dir, "05_report_data.RData")

# Build the CAR-ready table: everything needed for CAR except %OC / LOI
# Sedimentation rate comes from CRS; bulk density from the raw activity profile.
car_ready <- NULL
if (!is.null(crs_ages) && "sed_rate_cm_yr" %in% names(crs_ages) &&
    !is.null(activity_df)) {

  bd_df <- activity_df %>%
    mutate(depth_cm = (upper_depth_cm + lower_depth_cm) / 2) %>%
    select(depth_cm, bulk_density_g_cm3 = dry_bulk_density_g_cm3)

  car_ready <- crs_ages %>%
    filter(!is.na(age), !is.na(sed_rate_cm_yr), sed_rate_cm_yr > 0) %>%
    select(depth_cm, age_yr_BP = age, age_95_lo = age_min_yr,
           age_95_hi = age_max_yr, sed_rate_cm_yr) %>%
    left_join(bd_df, by = "depth_cm") %>%
    mutate(
      # pre-compute the non-carbon terms so only LOI/OC needs to be added later
      # CAR (g C m-2 yr-1) = C_fraction * BD (g/cm3) * SAR (cm/yr) * 10000
      car_factor = bulk_density_g_cm3 * sed_rate_cm_yr * 10000,
      ci_width_yr = age_95_hi - age_95_lo,
      ci_pct_age  = round(ci_width_yr / pmax(age_yr_BP, 1) * 100, 1),
      uncertainty_flag = ifelse(ci_pct_age > 50 | is.na(ci_pct_age),
                                "HIGH — more data needed", "acceptable")
    )

  write.csv(car_ready,
            file.path(output_dir, "05_car_ready_table.csv"),
            row.names = FALSE)
  cat("Saved: 05_car_ready_table.csv  (CAR-ready age-depth table)\n")
}

# Summarise Cs-137 availability for the report
has_cs137 <- !is.null(activity_df) && "cs137_dpm_g" %in% names(activity_df) &&
             any(!is.na(activity_df$cs137_dpm_g))

# Compute per-model CI summary for the report
ci_summary <- all_ages %>%
  mutate(ci_width = age_max_yr - age_min_yr,
         ci_pct   = ci_width / pmax(age, 1) * 100) %>%
  group_by(model) %>%
  summarise(
    n_sections       = n(),
    datable_depth_cm = max(depth_cm, na.rm = TRUE),
    oldest_age_yr    = max(age,      na.rm = TRUE),
    mean_ci_pct      = round(mean(ci_pct, na.rm = TRUE), 1),
    n_high_unc       = sum(ci_pct > 50, na.rm = TRUE),
    .groups = "drop"
  )

save(all_ages, activity_df, age_wide, ci_summary,
     crs_stats, rplum_stats, serac_stats, bayes_stats,
     car_ready, has_cs137, palette_colours, output_dir,
     file = report_data_path)
cat("Saved report data snapshot:", report_data_path, "\n")


# --- 9b. Write the .qmd file ---
qmd_path <- file.path(output_dir, "05_RERCA_Report.qmd")

qmd_text <- r"(---
title: "Pb-210 Sediment Chronology — RERCA Summary Report"
subtitle: "Multi-model age-depth comparison and carbon accumulation readiness assessment"
date: today
format:
  html:
    toc: true
    toc-depth: 3
    number-sections: true
    theme: cosmo
    embed-resources: true
    code-fold: true
    fig-width: 8
    fig-height: 5
execute:
  echo: false
  warning: false
  message: false
---

```{r setup}
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)
library(scales)

# Load snapshot written by Script 05
load("05_report_data.RData")

# Helper: flag colour
flag_colour <- function(x) {
  ifelse(x == "acceptable",
         '<span style="color:#1a9641;font-weight:bold;">✔ acceptable</span>',
         '<span style="color:#d7191c;font-weight:bold;">⚠ HIGH — more data needed</span>')
}
```

## Executive Summary

This report summarises the multi-model Pb-210 sediment chronology produced by
the RERCA workflow (scripts 01–05). Four independent age-depth methods were
compared: **pb210 CRS**, **rplum**, **serac CRS**, **serac CFCS**, and
**Bayesian (Stan)**. The report evaluates model agreement, flags sections where
chronological uncertainty is too large for reliable carbon accumulation
reconstruction, and prepares the age-depth dataset for carbon accumulation rate
(RERCA) calculation once organic carbon or LOI measurements are available.

> **Carbon data status:** No LOI or %OC measurements have been provided yet.
> A CAR-ready table (all non-carbon terms pre-calculated) has been saved as
> `05_car_ready_table.csv`. Add a `loi_percent` or `oc_percent` column and
> multiply by the `car_factor` column to compute RERCA.

---

## Pb-210 Activity Profile

```{r activity-plot, fig.cap="Total and unsupported Pb-210 activity with depth. The dashed line marks the estimated background (supported) activity level."}
if (!is.null(activity_df)) {
  bg_line <- unique(activity_df$pb210_supported)[1]

  act_long <- activity_df %>%
    select(depth_cm, total = pb210_total_dpm_g, unsupported = pb210_unsupported) %>%
    pivot_longer(c(total, unsupported), names_to = "fraction", values_to = "activity") %>%
    left_join(
      activity_df %>% select(depth_cm,
                             total_err = pb210_error_dpm_g,
                             unsup_err = pb210_unsupported_err),
      by = "depth_cm"
    ) %>%
    mutate(
      err      = ifelse(fraction == "total", total_err, unsup_err),
      fraction = recode(fraction,
                        total       = "Total Pb-210",
                        unsupported = "Unsupported Pb-210")
    )

  surface_act   <- activity_df$pb210_total_dpm_g[which.min(activity_df$depth_cm)]
  dynamic_range <- round(surface_act / bg_line, 1)

  ggplot(act_long, aes(x = activity, y = depth_cm, colour = fraction)) +
    geom_point(size = 2.5) +
    geom_errorbarh(aes(xmin = activity - err, xmax = activity + err), height = 0.3) +
    geom_vline(xintercept = bg_line, linetype = "dashed", colour = "grey40") +
    annotate("text", x = bg_line, y = max(activity_df$depth_cm) * 0.9,
             label = sprintf("Background\n%.3f DPM/g", bg_line),
             hjust = -0.1, size = 3.5, colour = "grey40") +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "Activity (DPM / g dry wt)") +
    scale_colour_manual(values = c("Total Pb-210"       = "#2c7bb6",
                                   "Unsupported Pb-210" = "#d7191c")) +
    labs(colour = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
} else {
  cat("Activity profile data not available.")
}
```

```{r activity-table}
if (!is.null(activity_df)) {
  surface_act   <- activity_df$pb210_total_dpm_g[which.min(activity_df$depth_cm)]
  bg_line       <- unique(activity_df$pb210_supported)[1]
  dynamic_range <- round(surface_act / bg_line, 1)

  dr_note <- dplyr::case_when(
    dynamic_range >= 5  ~ "Robust signal — strong model constraints expected.",
    dynamic_range >= 3  ~ "Moderate signal — models should be reliable but consider additional measurements.",
    TRUE                ~ "Weak signal — dating uncertainty will be high; additional Pb-210 measurements strongly recommended."
  )

  knitr::kable(data.frame(
    Metric  = c("Measured sections", "Background (supported)", "Surface total activity",
                "Dynamic range (surface/background)", "Assessment"),
    Value   = c(nrow(activity_df),
                sprintf("%.3f DPM/g", bg_line),
                sprintf("%.3f DPM/g", surface_act),
                sprintf("%.1fx", dynamic_range),
                dr_note)
  ), col.names = c("Metric", "Value"))
}
```

---

## Age-Depth Model Comparison

```{r age-depth-plot, fig.cap="All Pb-210 age-depth models overlaid. Shaded bands show 95% uncertainty intervals. Model selection is left to the analyst based on site knowledge and the diagnostics below."}
ggplot(all_ages, aes(x = age, y = depth_cm, colour = model, fill = model)) +
  geom_ribbon(aes(xmin = age_min_yr, xmax = age_max_yr),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Age (years before coring)") +
  scale_colour_manual(values = palette_colours) +
  scale_fill_manual(values = palette_colours, guide = "none") +
  labs(subtitle = "Shaded bands = 95% uncertainty intervals", colour = "Model") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))
```

### Age comparison at 2-cm depth intervals

```{r age-comparison-table}
age_wide %>%
  select(-starts_with("model_age_")) %>%
  mutate(across(where(is.numeric), \(x) round(x, 1))) %>%
  knitr::kable(caption = "Interpolated ages (yr before coring) at common 2-cm depth intervals.")
```

### Inter-model agreement

```{r agreement-table}
age_wide %>%
  select(depth_cm, model_age_min, model_age_max, model_age_range) %>%
  mutate(across(where(is.numeric), \(x) round(x, 1)),
         agreement = ifelse(model_age_range < 10, "good", "poor")) %>%
  knitr::kable(
    col.names = c("Depth (cm)", "Min age (yr)", "Max age (yr)",
                  "Spread (yr)", "Agreement"),
    caption = "Spread < 10 yr = good model agreement; > 10 yr warrants investigation."
  )
```

---

## Model Diagnostics

```{r diagnostics-table}
ci_summary %>%
  mutate(
    flag = ifelse(mean_ci_pct > 50, "⚠ wide", "✔ ok"),
    n_high_unc = paste0(n_high_unc, " / ", n_sections, " sections")
  ) %>%
  knitr::kable(
    col.names = c("Model", "Sections dated", "Datable depth (cm)",
                  "Oldest age (yr BP)", "Mean 95% CI (% of age)",
                  "Sections with CI > 50%", "CI flag"),
    caption = paste0("95% CI > 50% of estimated age flags sections where uncertainty",
                     " is too large for reliable RERCA. Additional measurements are",
                     " recommended at those depths.")
  )
```

```{r serac-cfcs}
if (!is.null(serac_stats) && "cfcs_r2" %in% names(serac_stats)) {
  cfcs_row <- serac_stats %>% filter(!is.na(cfcs_r2))
  if (nrow(cfcs_row) > 0) {
    r2    <- cfcs_row$cfcs_r2[1]
    sar   <- abs(cfcs_row$cfcs_sar_mm_yr[1])
    sar_e <- abs(cfcs_row$cfcs_sar_err[1])
    note  <- if (r2 > 0.9) {
      "R² > 0.9: log-linear decay is consistent with constant sedimentation. CFCS assumptions reasonable."
    } else if (r2 > 0.8) {
      "R² 0.8–0.9: moderate fit. CRS is likely more appropriate than CFCS."
    } else {
      "R² < 0.8: poor fit. Sedimentation rate is variable — CRS or Bayesian model preferred."
    }
    knitr::kable(data.frame(
      Metric = c("CFCS SAR (mm/yr)", "CFCS R²", "Interpretation"),
      Value  = c(sprintf("%.3f ± %.3f", sar, sar_e), round(r2, 3), note)
    ), col.names = c("Metric", "Value"),
    caption = "CFCS linearity check (serac). R² measures how well a single constant sedimentation rate explains the Pb-210 decay profile.")
  }
}
```

```{r bayesian-diagnostics}
if (!is.null(bayes_stats) && "n_divergent" %in% names(bayes_stats)) {
  divs <- bayes_stats$n_divergent[1]
  div_note <- if (divs == 0) {
    "0 divergences — sampler converged well."
  } else if (divs <= 10) {
    sprintf("%d divergences — minor; consider tightening priors.", divs)
  } else {
    sprintf("%d divergences — poor MCMC mixing. Results should be treated cautiously.", divs)
  }
  knitr::kable(data.frame(
    Metric = c("Divergent transitions", "phi (supported) mean",
               "phi 95% CI", "Assessment"),
    Value  = c(divs,
               sprintf("%.3f DPM/g", bayes_stats$phi_mean_dpm_g[1]),
               sprintf("%.3f – %.3f DPM/g", bayes_stats$phi_95_lo[1], bayes_stats$phi_95_hi[1]),
               div_note)
  ), col.names = c("Metric", "Value"),
  caption = "Bayesian (Stan) MCMC diagnostics.")
}
```

---

## Cs-137 Validation

```{r cs137-section, results='asis'}
if (!has_cs137) {
  cat("> **Cs-137 data not present in this dataset.**\n>\n",
      "> Adding Cs-137 measurements is strongly recommended as an independent\n",
      "> time marker. The 1963 peak (peak atmospheric fallout) provides a fixed\n",
      "> date to validate or anchor the Pb-210 age-depth models.\n>\n",
      "> **What to do:** Include `cs137_dpm_g` and `cs137_error_dpm_g` columns\n",
      "> in your input CSV and re-run scripts 01–05.\n")
}
```

```{r cs137-plot, eval=has_cs137, fig.cap="Cs-137 activity profile. The 1963 peak (maximum atmospheric fallout) provides an independent time marker for validating the Pb-210 age-depth models."}
if (has_cs137) {
  cs_df <- activity_df %>%
    filter(!is.na(cs137_dpm_g)) %>%
    mutate(depth_cm = (upper_depth_cm + lower_depth_cm) / 2)

  peak_depth <- cs_df$depth_cm[which.max(cs_df$cs137_dpm_g)]

  # Look up what age each model assigns to the peak depth
  age_at_peak <- all_ages %>%
    group_by(model) %>%
    summarise(
      age_at_cs137 = approx(depth_cm, age, xout = peak_depth, rule = 1)$y,
      .groups = "drop"
    ) %>%
    filter(!is.na(age_at_cs137)) %>%
    mutate(
      year_AD = lubridate::year(Sys.Date()) - age_at_cs137,
      expected_yr = 1963,
      offset_yr   = round(year_AD - expected_yr, 1)
    )

  p_cs <- ggplot(cs_df, aes(x = cs137_dpm_g, y = depth_cm)) +
    geom_point(size = 2.5, colour = "#756bb1") +
    geom_errorbarh(aes(xmin = cs137_dpm_g - cs137_error_dpm_g,
                       xmax = cs137_dpm_g + cs137_error_dpm_g), height = 0.3,
                   colour = "#756bb1") +
    geom_hline(yintercept = peak_depth, linetype = "dashed", colour = "grey40") +
    annotate("text", x = max(cs_df$cs137_dpm_g) * 0.7, y = peak_depth,
             label = sprintf("Peak at %.1f cm", peak_depth),
             vjust = -0.5, size = 3.5, colour = "grey40") +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "Cs-137 activity (DPM / g dry wt)") +
    theme_bw(base_size = 13)
  print(p_cs)

  knitr::kable(age_at_peak,
    col.names = c("Model", "Age at Cs-137 peak (yr BP)",
                  "Modelled year (AD)", "Expected year (AD)", "Offset (yr)"),
    caption = sprintf(
      "Cs-137 peak detected at %.1f cm depth. Each model's assigned age is compared to the expected 1963 date. Offset < 5 yr = good validation; > 10 yr warrants investigation.", peak_depth))
}
```

---

## Carbon Accumulation Rate (RERCA) Readiness

### What is RERCA?

Recent Rates of Carbon Accumulation (RERCA) are calculated as:

$$
\text{RERCA} \; (\text{g C m}^{-2} \text{yr}^{-1}) = C_{\text{fraction}} \times \rho_b \; (\text{g/cm}^3) \times \text{SAR} \; (\text{cm/yr}) \times 10{,}000
$$

where $C_{\text{fraction}}$ is the proportion of organic carbon (from LOI or elemental analysis), $\rho_b$ is dry bulk density, and SAR is the sedimentation accumulation rate from the age-depth model.

### CAR-ready table (carbon data pending)

All terms except the carbon fraction have been pre-calculated. The `car_factor`
column equals $\rho_b \times \text{SAR} \times 10{,}000$ — multiply it by your
LOI-derived or measured carbon fraction to get RERCA.

```{r car-table}
if (!is.null(car_ready)) {
  car_ready %>%
    mutate(
      age_yr_BP      = round(age_yr_BP, 1),
      age_95_lo      = round(age_95_lo, 1),
      age_95_hi      = round(age_95_hi, 1),
      sed_rate_cm_yr = round(sed_rate_cm_yr, 4),
      bulk_density_g_cm3 = round(bulk_density_g_cm3, 3),
      car_factor     = round(car_factor, 2),
      ci_pct_age     = paste0(ci_pct_age, "%")
    ) %>%
    knitr::kable(
      col.names = c("Depth (cm)", "Age (yr BP)", "95% lo", "95% hi",
                    "SAR (cm/yr)", "Bulk density (g/cm³)",
                    "CAR factor", "CI (% of age)", "Uncertainty flag"),
      caption = paste0(
        "CAR-ready table from CRS model. car_factor = bulk_density × SAR × 10000.",
        " Multiply by C_fraction (LOI/OC) to get RERCA (g C m⁻² yr⁻¹).",
        " Flagged rows have 95% CI > 50% of age and require additional measurements.")
    )
} else {
  cat("> CAR-ready table could not be built — CRS ages or activity profile missing.\n")
}
```

```{r car-flag-summary, results='asis'}
if (!is.null(car_ready)) {
  n_high <- sum(car_ready$uncertainty_flag != "acceptable", na.rm = TRUE)
  n_ok   <- sum(car_ready$uncertainty_flag == "acceptable", na.rm = TRUE)
  cat(sprintf(
    "\n**%d of %d sections** have acceptable uncertainty (95%% CI < 50%% of age).",
    n_ok, nrow(car_ready)))
  if (n_high > 0) {
    bad_depths <- car_ready$depth_cm[car_ready$uncertainty_flag != "acceptable"]
    cat(sprintf(
      " **%d section(s) flagged** at depths: %s cm — additional measurements needed before these can be used for RERCA.\n",
      n_high, paste(round(bad_depths, 1), collapse = ", ")))
  }
}
```

---

## What More Is Needed

The following measurements are required or recommended to complete a reliable
RERCA reconstruction. Items are ordered by priority.

### Required — not yet available

| Priority | Measurement | Why needed | Column in template |
|----------|-------------|------------|-------------------|
| **1** | LOI (%) or %OC (elemental) | The carbon fraction in the RERCA equation — without this, CAR cannot be calculated | `loi_percent` / add `oc_percent` |
| **2** | Carbon content for *all* dated sections | RERCA requires continuous coverage — gaps create unresolvable bias | Same as above |

### Strongly recommended — to reduce chronological uncertainty

| Priority | Measurement | Why needed | Column in template |
|----------|-------------|------------|-------------------|
| **3** | Additional Pb-210 measurements at depths with high uncertainty | Reduces CI width in flagged sections (see table above) | `pb210_total_dpm_g` / `pb210_error_dpm_g` |
| **4** | Cs-137 activity profile (if not already measured) | Independent 1963 time marker to validate or anchor the Pb-210 model | `cs137_dpm_g` / `cs137_error_dpm_g` |
| **5** | Ra-226 measurements at multiple depths | Direct measurement of supported (background) Pb-210; currently estimated from deep-core tail — Ra-226 reduces this uncertainty significantly | `ra226_supported_dpm_g` / `ra226_error_dpm_g` |

### Optional — to strengthen the reconstruction

| Priority | Measurement | Why needed |
|----------|-------------|------------|
| **6** | Spheroidal carbonaceous particles (SCPs) | Additional stratigraphic markers (~1850, ~1950, ~1980 peaks) for independent age constraints |
| **7** | Replicate Pb-210 measurements on select sections | Quantify analytical reproducibility, especially in sections near the detection limit |
| **8** | Core top and bottom Pb-210 if not measured | Ensures the model is anchored at the sediment-water interface |

---

## Next Steps

1. **Select a preferred age-depth model** from the comparison plot and diagnostic
   tables above. Consider CRS as the default for variable-sedimentation sites;
   use CFCS only if R² > 0.9 (constant sedimentation confirmed).

2. **Collect organic carbon data** (LOI% or %OC by elemental analysis) for all
   sections listed in the CAR-ready table. If using LOI, apply a site-appropriate
   conversion factor (e.g. Craft et al. 1991 for coastal wetlands; Heiri et al.
   2001 for lake sediments).

3. **Address flagged sections** — obtain additional Pb-210 measurements or
   Cs-137 / Ra-226 data at depths where uncertainty is > 50% of the estimated
   age before including those sections in the RERCA calculation.

4. **Add Cs-137 data** and re-run scripts 01–05. The validation section of this
   report will automatically populate once `cs137_dpm_g` values are present in
   the input CSV.

5. **Calculate RERCA** by loading `05_car_ready_table.csv`, joining your carbon
   data, and computing: `RERCA = C_fraction * car_factor`.

6. **Consider uncertainty propagation** — propagate the 95% CI on age through
   to the RERCA estimates, particularly for sections flagged as high-uncertainty.

---

## Session Info

```{r session-info}
sessionInfo()
```
)"

writeLines(qmd_text, qmd_path)
cat("Saved Quarto report template:", qmd_path, "\n")


# --- 9c. Render the report ---
cat("Rendering Quarto report...\n")
tryCatch({
  quarto::quarto_render(
    input  = qmd_path,
    output_format = "html"
  )
  html_out <- sub("\\.qmd$", ".html", qmd_path)
  cat(sprintf("Report rendered: %s\n", html_out))
}, error = function(e) {
  cat("Could not render automatically:", conditionMessage(e), "\n")
  cat("To render manually, run:\n")
  cat(sprintf('  quarto::quarto_render("%s")\n', qmd_path))
  cat("  or from the terminal: quarto render", qmd_path, "\n")
})
