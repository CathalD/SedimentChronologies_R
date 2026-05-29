# =============================================================================
# RERCA | Script 01: Pb-210 CRS Model using the {pb210} package
# Method: Constant Rate of Supply (CRS)
# Package: https://github.com/paleolimbot/pb210
# =============================================================================
#
# OVERVIEW
# --------
# The CRS model assumes Pb-210 is supplied to the sediment surface at a
# constant rate over time, but allows sedimentation rate to vary. It is the
# most widely used model for lake sediment Pb-210 chronologies.
#
# This script requires only measured TOTAL Pb-210 activity. Supported
# (background) activity is estimated from the data if not directly measured.
#
# STEPS
# -----
# 1. Install / load packages
# 2. Load your data
# 3. Interpolate missing bulk density values
# 4. Estimate supported (background) Pb-210
# 5. Calculate unsupported (excess) Pb-210
# 6. Fit the CRS model
# 7. Plot results
# 8. Export age-depth table to console
# =============================================================================


# =============================================================================
# STEP 1 — Install and load packages
# =============================================================================
# Run this block once. After packages are installed, you can comment out
# the install lines.

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("pb210",   quietly = TRUE)) remotes::install_github("paleolimbot/pb210")
if (!requireNamespace("errors",  quietly = TRUE)) install.packages("errors")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")

library(pb210)
library(errors)
library(dplyr)
library(ggplot2)
library(tidyr)


# =============================================================================
# STEP 2 — Load your data
# =============================================================================
# Point data_file at your filled-in CSV (see data/template_pb210_data.csv).
# The example dataset ships with this repository.

data_file <- "RERCA/data/example_pb210_data.csv"   # <-- CHANGE THIS: paste the full path to your CSV here, e.g.:

core_raw <- read.csv(data_file, na.strings = c("NA", ""))

# Quick sanity check — print the first few rows
cat("\n--- Raw data preview ---\n")
print(head(core_raw, 10))
cat(sprintf("\nRows loaded: %d\n", nrow(core_raw)))

# Required columns
required_cols <- c("upper_depth_cm", "lower_depth_cm",
                   "dry_bulk_density_g_cm3",
                   "pb210_total_dpm_g", "pb210_error_dpm_g")
missing_cols <- setdiff(required_cols, names(core_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
       "\nCheck your CSV matches the template column names.")
}

# Add midpoint depth column
core_raw <- core_raw %>%
  mutate(depth_cm = (upper_depth_cm + lower_depth_cm) / 2)


# =============================================================================
# STEP 3 — Interpolate missing bulk density values
# =============================================================================
# Many labs only measure bulk density on a subset of sections. Missing values
# are linearly interpolated here so that cumulative mass depth can be
# calculated for all intervals.

if (any(is.na(core_raw$dry_bulk_density_g_cm3))) {
  cat("\nNote: interpolating", sum(is.na(core_raw$dry_bulk_density_g_cm3)),
      "missing bulk density value(s) by linear interpolation.\n")
  core_raw <- core_raw %>%
    mutate(dry_bulk_density_g_cm3 = approx(
      x    = depth_cm[!is.na(dry_bulk_density_g_cm3)],
      y    = dry_bulk_density_g_cm3[!is.na(dry_bulk_density_g_cm3)],
      xout = depth_cm,
      rule = 2   # extrapolate at edges using nearest value
    )$y)
}

# Calculate slice thickness and cumulative mass depth (g/cm²)
core_raw <- core_raw %>%
  mutate(
    slice_thickness_cm  = lower_depth_cm - upper_depth_cm,
    mass_depth_g_cm2    = cumsum(dry_bulk_density_g_cm3 * slice_thickness_cm)
  )

cat("\n--- Processed data (depths & mass depth) ---\n")
print(core_raw %>% select(upper_depth_cm, lower_depth_cm, depth_cm,
                           dry_bulk_density_g_cm3, mass_depth_g_cm2))


# =============================================================================
# STEP 4 — Estimate supported (background) Pb-210
# =============================================================================
# Supported Pb-210 is the activity contributed by in-situ decay of Ra-226 in
# the sediment matrix. It must be subtracted to get the unsupported (excess)
# fraction that carries the age signal.
#
# TWO OPTIONS:
#   A) Use measured Ra-226 values (if available in your dataset)
#   B) Estimate background from the asymptotic tail of the total activity
#      profile (appropriate when Ra-226 was not measured)
#
# The script checks for Option A first and falls back to Option B if Ra-226
# data are absent or mostly missing.

has_ra226 <- "ra226_supported_dpm_g" %in% names(core_raw) &&
             sum(!is.na(core_raw$ra226_supported_dpm_g)) >= 3

if (has_ra226) {
  # --- Option A: use measured Ra-226 ---
  cat("\nSupported Pb-210 source: measured Ra-226 values.\n")
  core_raw <- core_raw %>%
    mutate(
      pb210_supported      = ra226_supported_dpm_g,
      pb210_supported_err  = ifelse("ra226_error_dpm_g" %in% names(.),
                                    ra226_error_dpm_g, 0)
    )
} else {
  # --- Option B: asymptotic tail estimation ---
  cat("\nNo Ra-226 measurements found (or fewer than 3 values).\n")
  cat("Estimating supported Pb-210 from the deep-section tail.\n\n")

  # Ask the user whether to use the automatic asymptotic estimate or supply
  # a manual value.
  cat("Choose background estimation method:\n")
  cat("  1 = Automatic: use the mean of the lowest-activity sections\n")
  cat("  2 = Manual: I will enter my own background value\n")
  method_choice <- readline(prompt = "Enter 1 or 2: ")

  if (trimws(method_choice) == "2") {
    bg_value <- as.numeric(readline(
      prompt = "Enter supported Pb-210 activity (DPM/g dry wt): "))
    bg_error <- as.numeric(readline(
      prompt = "Enter uncertainty on that value (DPM/g, 1 SD): "))
    cat(sprintf("\nUsing user-supplied background: %.3f ± %.3f DPM/g\n",
                bg_value, bg_error))
  } else {
    # Automatic: take mean of the bottom 20% of sections that have activity data
    measured_rows <- core_raw %>% filter(!is.na(pb210_total_dpm_g))
    n_tail <- max(2, round(nrow(measured_rows) * 0.2))
    tail_rows <- tail(measured_rows, n_tail)
    bg_value <- mean(tail_rows$pb210_total_dpm_g, na.rm = TRUE)
    bg_error <- sd(tail_rows$pb210_total_dpm_g,   na.rm = TRUE)
    if (is.na(bg_error) || bg_error == 0) bg_error <- bg_value * 0.1  # fallback 10%
    cat(sprintf(
      "\nAutomatic background estimate: %.3f ± %.3f DPM/g\n(mean of %d deepest measured sections)\n",
      bg_value, bg_error, n_tail))
    cat("If this looks unreasonable, re-run and choose option 2 to enter manually.\n")
  }

  core_raw <- core_raw %>%
    mutate(
      pb210_supported     = bg_value,
      pb210_supported_err = bg_error
    )
}


# =============================================================================
# STEP 5 — Calculate unsupported (excess) Pb-210
# =============================================================================

core_raw <- core_raw %>%
  mutate(
    pb210_unsupported     = pb210_total_dpm_g  - pb210_supported,
    pb210_unsupported_err = sqrt(pb210_error_dpm_g^2 + pb210_supported_err^2)
  )

# Flag and report sections where unsupported activity is not significantly
# above background (these are near or below the detection limit)
core_raw <- core_raw %>%
  mutate(above_background = pb210_unsupported > pb210_unsupported_err * 1.5)

n_below <- sum(!core_raw$above_background & !is.na(core_raw$pb210_unsupported))
if (n_below > 0) {
  cat(sprintf(
    "\nNote: %d section(s) have unsupported Pb-210 not significantly above background.\n",
    n_below))
  cat("These sections are below the datable horizon and will be excluded from the CRS model.\n")
}

cat("\n--- Unsupported Pb-210 profile ---\n")
print(core_raw %>%
  filter(!is.na(pb210_total_dpm_g)) %>%
  select(depth_cm, pb210_total_dpm_g, pb210_supported, pb210_unsupported,
         pb210_unsupported_err, above_background))


# =============================================================================
# STEP 6 — Fit the CRS model using the pb210 package
# =============================================================================
# Workflow (following the pb210 package vignette):
#   1. pb210_excess()          — subtract background from total activity
#   2. pb210_cumulative_mass() — running sum of mass per unit area (g/cm²)
#   3. pb210_inventory()       — total unsupported Pb-210 inventory
#   4. pb210_crs() + predict() — CRS age model

cat("\n--- Fitting CRS model ---\n")

# Step 6a: pack total activity and background into errors-class vectors,
# then compute unsupported (excess) activity via pb210_excess()
total_activity_with_err <- set_errors(
  core_raw$pb210_total_dpm_g,
  core_raw$pb210_error_dpm_g
)
background_with_err <- set_errors(bg_value, bg_error)

core_raw$excess_pb210 <- pb210_excess(total_activity_with_err, background_with_err)

# Step 6b: cumulative dry mass per unit area (g/cm²) from bulk density × thickness
core_raw$cumulative_dry_mass <- pb210_cumulative_mass(
  core_raw$dry_bulk_density_g_cm3 * core_raw$slice_thickness_cm
)

# Step 6c: calculate total inventory — fitted over all rows with measured activity
measured <- core_raw %>% filter(!is.na(pb210_total_dpm_g))

inventory <- pb210_inventory(
  measured$cumulative_dry_mass,
  measured$excess_pb210,
  model_bottom = ~pb210_fit_loglinear(
    ..1, ..2,
    subset = ~finite_tail(..1, ..2, n_tail = 2)
  )
)

# Step 6d: fit CRS model and extract predicted ages
crs_ages <- pb210_crs(
  measured$cumulative_dry_mass,
  measured$excess_pb210,
  inventory = inventory
) %>%
  predict()

cat("\nRaw CRS output columns:", paste(names(crs_ages), collapse = ", "), "\n")
print(crs_ages)

# Attach ages to the measured sections
core_dated <- measured %>%
  bind_cols(crs_ages) %>%
  mutate(
    age_min_yr = age - age_sd * 1.96,
    age_max_yr = age + age_sd * 1.96
  )

# Print age-depth table
cat("\n=== CRS Age-Depth Results ===\n")
print(core_dated %>%
  select(depth_cm, age, age_sd, age_min_yr, age_max_yr) %>%
  mutate(across(where(is.numeric), \(x) round(x, 1))))

cat(sprintf(
  "\nDatable horizon: %.1f cm  |  oldest age: %.0f yr BP (95%% CI: %.0f–%.0f yr)\n",
  max(core_dated$depth_cm),
  max(core_dated$age,     na.rm = TRUE),
  max(core_dated$age_min_yr, na.rm = TRUE),
  max(core_dated$age_max_yr, na.rm = TRUE)
))


# =============================================================================
# STEP 7 — Plots
# =============================================================================

# --- Plot A: Total and unsupported Pb-210 activity profile ---
pb210_profile_data <- core_raw %>%
  filter(!is.na(pb210_total_dpm_g)) %>%
  select(depth_cm, pb210_total_dpm_g, pb210_error_dpm_g,
         pb210_unsupported, pb210_unsupported_err) %>%
  pivot_longer(cols = c(pb210_total_dpm_g, pb210_unsupported),
               names_to  = "fraction",
               values_to = "activity") %>%
  mutate(
    error = ifelse(fraction == "pb210_total_dpm_g",
                   pb210_error_dpm_g, pb210_unsupported_err),
    fraction = recode(fraction,
      pb210_total_dpm_g    = "Total Pb-210",
      pb210_unsupported    = "Unsupported Pb-210"
    )
  )

p1 <- ggplot(pb210_profile_data,
             aes(x = activity, y = depth_cm, colour = fraction)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = activity - error, xmax = activity + error),
                 height = 0.3) +
  geom_vline(xintercept = unique(core_raw$pb210_supported[!is.na(core_raw$pb210_supported)]),
             linetype = "dashed", colour = "grey50") +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Activity (DPM / g dry wt)") +
  scale_colour_manual(values = c("Total Pb-210" = "#2c7bb6",
                                 "Unsupported Pb-210" = "#d7191c")) +
  labs(title = "Pb-210 Activity Profile",
       subtitle = "Dashed line = estimated supported (background) activity",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p1)

# --- Plot B: CRS age-depth model ---
p2 <- ggplot(core_dated, aes(x = age, y = depth_cm)) +
  geom_ribbon(aes(xmin = age_min_yr, xmax = age_max_yr),
              fill = "#2c7bb6", alpha = 0.25) +
  geom_line(colour = "#2c7bb6", linewidth = 1) +
  geom_point(colour = "#2c7bb6", size = 2.5) +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Age (years before coring)") +
  labs(title   = "CRS Age-Depth Model",
       subtitle = "Shaded band = 95% confidence interval") +
  theme_bw(base_size = 12)

print(p2)

# --- Plot C: Sedimentation rate ---
core_dated <- core_dated %>%
  arrange(depth_cm) %>%
  mutate(sed_rate_cm_yr = c(NA, diff(depth_cm) / diff(age)))

p3 <- ggplot(core_dated %>% filter(!is.na(sed_rate_cm_yr)),
             aes(x = sed_rate_cm_yr, y = depth_cm)) +
  geom_step(orientation = "y", colour = "#1a9641", linewidth = 1) +
  geom_point(colour = "#1a9641", size = 2.5) +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Sedimentation rate (cm / yr)") +
  labs(title = "CRS Sedimentation Rate Profile") +
  theme_bw(base_size = 12)

print(p3)

cat("\nScript 01 complete. Review plots and the age-depth table above.\n")
cat("Next: run 02_rplum.R for the Bayesian Pb-210 approach.\n")
