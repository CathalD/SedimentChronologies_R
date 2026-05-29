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
if (!is.null(crs_ages) && "sed_rate_cm_yr" %in% names(crs_ages)) {

  car_ready <- crs_ages %>%
    filter(!is.na(age), !is.na(sed_rate_cm_yr), sed_rate_cm_yr > 0) %>%
    select(depth_cm, age_yr_BP = age, age_95_lo = age_min_yr,
           age_95_hi = age_max_yr, sed_rate_cm_yr)

  # Join bulk density from the activity profile if the column was saved there
  if (!is.null(activity_df) && "dry_bulk_density_g_cm3" %in% names(activity_df)) {
    bd_df <- activity_df %>%
      select(depth_cm, bulk_density_g_cm3 = dry_bulk_density_g_cm3)
    car_ready <- car_ready %>% left_join(bd_df, by = "depth_cm")
  } else {
    car_ready <- car_ready %>% mutate(bulk_density_g_cm3 = NA_real_)
    cat("Note: bulk density not found in activity profile CSV.\n")
    cat("      Add a bulk_density_g_cm3 column to 05_car_ready_table.csv manually.\n")
  }

  car_ready <- car_ready %>%
    mutate(
      # CAR (g C m-2 yr-1) = C_fraction * BD (g/cm3) * SAR (cm/yr) * 10000
      car_factor       = bulk_density_g_cm3 * sed_rate_cm_yr * 10000,
      ci_width_yr      = age_95_hi - age_95_lo,
      ci_pct_age       = round(ci_width_yr / pmax(age_yr_BP, 1) * 100, 1),
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


# --- 9b. Copy the .qmd to the output directory and render ---
# The .qmd is a standalone file (RERCA/05_RERCA_Report.qmd) so there are no
# quoting issues. It loads 05_report_data.RData from the same directory.
qmd_source <- file.path(dirname(output_dir), "05_RERCA_Report.qmd")
qmd_dest   <- file.path(output_dir, "05_RERCA_Report.qmd")

if (file.exists(qmd_source)) {
  file.copy(qmd_source, qmd_dest, overwrite = TRUE)
  cat("Copied report template to:", qmd_dest, "\n")
} else {
  cat("Warning: 05_RERCA_Report.qmd not found at", qmd_source, "\n")
  cat("Make sure it is in the RERCA/ directory alongside the other scripts.\n")
}

# --- 9c. Render the report ---
if (file.exists(qmd_dest)) {
  cat("Rendering Quarto report...\n")
  tryCatch({
    quarto::quarto_render(
      input         = qmd_dest,
      output_format = "html"
    )
    html_out <- sub("\\.qmd$", ".html", qmd_dest)
    cat(sprintf("Report rendered: %s\n", html_out))
  }, error = function(e) {
    cat("Could not render automatically:", conditionMessage(e), "\n")
    cat("To render manually, open RStudio and click Render on:\n")
    cat(" ", qmd_dest, "\n")
    cat("Or run: quarto render", qmd_dest, "\n")
  })
}

