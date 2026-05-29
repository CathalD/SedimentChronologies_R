# =============================================================================
# RERCA | Script 02: Bayesian Pb-210 Chronology using {rplum}
# Method: Full Bayesian age-depth modelling with Plum
# Package: https://CRAN.R-project.org/package=rplum
# =============================================================================
#
# OVERVIEW
# --------
# rplum extends the Bacon age-depth modelling framework to incorporate
# Pb-210 data directly, jointly estimating the supported activity and
# sedimentation history within a Bayesian framework. It is especially
# powerful when you want to integrate Pb-210 with C-14 dates later.
#
# IMPORTANT: rplum writes output files to a directory called "Plum_runs/"
# in your working directory. Make sure your working directory is set to
# the project root (SedimentChronology_Workflow/).
#
# STEPS
# -----
# 1. Install / load packages
# 2. Prepare the rplum input file
# 3. Run the Plum model
# 4. Inspect and plot results
# 5. Extract the age-depth table
# =============================================================================


# =============================================================================
# STEP 1 — Install and load packages
# =============================================================================

if (!requireNamespace("rplum",   quietly = TRUE)) install.packages("rplum")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(rplum)
library(dplyr)
library(ggplot2)

# Set working directory to project root if not already there.
# Uncomment and adjust the line below if needed:
# setwd("path/to/SedimentChronology_Workflow")


# =============================================================================
# STEP 2 — Prepare the rplum input file
# =============================================================================
# rplum expects a specific CSV format with exactly these columns (in order):
#   depth, density, thickness, Pb210, sdPb210, (optional) Pb210supp, sdPb210supp
#
# The function below reads your standard project CSV and converts it.

data_file <- "RERCA/data/example_pb210_data.csv"   # <-- CHANGE THIS: paste the full path to your CSV here, e.g.:
# data_file <- "/Users/cathaldoherty/Desktop/CoastalBC_WorkflowV1/example_pb210_data.csv"
core_name <- "MyCore"                               # <-- CHANGE to a short core ID (no spaces)

core_raw <- read.csv(data_file, na.strings = c("NA", ""))

core_raw <- core_raw %>%
  mutate(
    depth     = (upper_depth_cm + lower_depth_cm) / 2,
    thickness = lower_depth_cm - upper_depth_cm,
    density   = dry_bulk_density_g_cm3
  )

# Interpolate missing bulk density
if (any(is.na(core_raw$density))) {
  core_raw <- core_raw %>%
    mutate(density = approx(
      x    = depth[!is.na(density)],
      y    = density[!is.na(density)],
      xout = depth, rule = 2)$y)
}

# Build the rplum input data frame
# rplum requires 6–9 columns in this exact order:
#   labID, depth, density, thickness, Pb210, sdPb210
#   (optional cols 7–9: Pb210supp, sdPb210supp, and one more optional field)
has_supported <- "ra226_supported_dpm_g" %in% names(core_raw) &&
                 sum(!is.na(core_raw$ra226_supported_dpm_g)) >= 3

if (has_supported) {
  plum_input <- core_raw %>%
    filter(!is.na(pb210_total_dpm_g)) %>%
    transmute(
      labID        = paste0(core_name, "_", seq_len(n())),
      depth        = depth,
      density      = density,
      thickness    = thickness,
      Pb210        = pb210_total_dpm_g,
      sdPb210      = pb210_error_dpm_g,
      Pb210supp    = ra226_supported_dpm_g,
      sdPb210supp  = ifelse("ra226_error_dpm_g" %in% names(.),
                            ra226_error_dpm_g, ra226_supported_dpm_g * 0.05)
    )
  cat("Supported Pb-210 column included from Ra-226 measurements.\n")
} else {
  plum_input <- core_raw %>%
    filter(!is.na(pb210_total_dpm_g)) %>%
    transmute(
      labID     = paste0(core_name, "_", seq_len(n())),
      depth     = depth,
      density   = density,
      thickness = thickness,
      Pb210     = pb210_total_dpm_g,
      sdPb210   = pb210_error_dpm_g
    )
  cat("No Ra-226 data found. rplum will estimate supported activity from the data.\n")
}

# Create the Plum_runs directory and write the input file
run_dir <- file.path("Plum_runs", core_name)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
plum_csv_path <- file.path(run_dir, paste0(core_name, ".csv"))
write.csv(plum_input, plum_csv_path, row.names = FALSE, quote = FALSE)

cat(sprintf("\nrplum input file written to: %s\n", plum_csv_path))
cat("Preview of rplum input:\n")
print(head(plum_input, 10))


# =============================================================================
# STEP 3 — Run the Plum model
# =============================================================================
# The Plum() function runs the MCMC sampler. This typically takes 1–5 minutes.
# Key parameters you may want to adjust:
#   thick     = section thickness for the age model (default 1 cm)
#   acc.mean  = prior mean for accumulation rate (yr/cm), adjust if needed
#   mem.mean  = prior for memory (autocorrelation between sections, 0–1)
#   burnin    = number of MCMC iterations to discard as burn-in

cat("\nRunning Plum MCMC sampler — this may take several minutes...\n")

Plum(
  core       = core_name,
  coredir    = "Plum_runs",
  thick      = 1,         # age-model section thickness in cm
  acc.mean   = 10,        # prior mean accumulation rate (yr/cm)
  acc.shape  = 1.5,       # prior shape for accumulation rate
  mem.mean   = 0.7,       # prior mean for memory
  mem.strength = 4,
  burnin     = 200,
  ssize      = 2000,      # number of posterior samples to keep
  BCAD       = FALSE,     # use years BP (before coring), set TRUE for AD/BC
  plot.pdf   = FALSE      # suppress automatic PDF (we plot manually below)
)

cat("\nPlum run complete.\n")


# =============================================================================
# STEP 4 — Inspect and plot results
# =============================================================================
# rplum creates diagnostic plots automatically. We also pull the age
# estimates manually for a clean ggplot.

# Plum stores results as <core_name>_<nsections>_ages.txt
# The section count varies, so we find the file by pattern match.
ages_file <- list.files(run_dir,
                        pattern = paste0("^", core_name, "_[0-9]+_ages\\.txt$"),
                        full.names = TRUE)

if (length(ages_file) == 0) {
  # Fallback: also try without section number (older rplum versions)
  ages_file <- file.path(run_dir, paste0(core_name, "_ages.txt"))
  if (!file.exists(ages_file)) ages_file <- character(0)
} else {
  ages_file <- ages_file[1]  # take the most recent if multiple exist
}

if (length(ages_file) > 0 && file.exists(ages_file)) {
  ages_df <- read.table(ages_file, header = TRUE)

  cat("\n=== rplum Posterior Age Estimates ===\n")
  print(ages_df %>% mutate(across(where(is.numeric), \(x) round(x, 1))))

  # Age-depth plot
  p_rplum <- ggplot(ages_df, aes(x = mean, y = depth)) +
    geom_ribbon(aes(xmin = min.95, xmax = max.95), fill = "#fdae61", alpha = 0.4) +
    geom_line(colour = "#d7191c", linewidth = 1) +
    geom_point(colour = "#d7191c", size = 2.5) +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "Age (years before coring)") +
    labs(title   = paste("rplum Bayesian Age-Depth Model:", core_name),
         subtitle = "Orange band = 95% CI; wider band = full posterior range") +
    theme_bw(base_size = 12)

  print(p_rplum)

} else {
  cat(sprintf("\nCould not find ages file in: %s\n", run_dir))
  cat("Files present:\n")
  print(list.files(run_dir))
  cat("Check that Plum() ran without errors above, and that core_name matches the directory.\n")
}

# Show the supported Pb-210 posterior estimate
cat("\nrplum estimates of supported Pb-210 activity are printed in the model summary above.\n")
cat("Look for 'phi' (the supported activity parameter) in the MCMC diagnostics.\n")


# =============================================================================
# STEP 5 — Compare with CRS (optional)
# =============================================================================
# If you ran 01_pb210_CRS.R, load its results and overlay them here.

crs_result_file <- "RERCA/output/crs_ages.csv"  # set this if you saved CRS output
if (file.exists(crs_result_file) && exists("ages_df")) {
  crs_ages <- read.csv(crs_result_file)
  p_compare <- ggplot() +
    geom_ribbon(data = ages_df,
                aes(x = mean, y = depth, xmin = X2.5., xmax = X97.5.),
                fill = "#fdae61", alpha = 0.4) +
    geom_line(data = ages_df,
              aes(x = mean, y = depth, colour = "rplum"), linewidth = 1) +
    geom_line(data = crs_ages,
              aes(x = age_yr, y = depth_cm, colour = "CRS"), linewidth = 1) +
    scale_y_reverse(name = "Depth (cm)") +
    scale_x_continuous(name = "Age (years before coring)") +
    scale_colour_manual(values = c("rplum" = "#d7191c", "CRS" = "#2c7bb6")) +
    labs(title = "CRS vs rplum Age-Depth Comparison", colour = "Model") +
    theme_bw(base_size = 12)
  print(p_compare)
}

cat("\nScript 02 complete.\n")
cat("Next: run 03_serac.R to apply CRS/CIC/CFCS models from the serac package.\n")
