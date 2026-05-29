# =============================================================================
# RERCA | Script 03: Pb-210 Chronology using {serac}
# Methods: CRS and CFCS models (compared side-by-side)
# Package: https://github.com/rosalieb/serac
# Reference: Bruel & Sabatier (2020), J. Environmental Radioactivity
# =============================================================================
#
# OVERVIEW
# --------
# serac (Short-livEd RAdioisotope Chronology) implements classical Pb-210
# dating models. This script runs CRS and CFCS on your core data and
# produces a side-by-side comparison of the two age-depth models.
#
#   CRS  — Constant Rate of Supply: assumes constant Pb-210 flux, allows
#           sedimentation rate to vary. Most widely used.
#
#   CFCS — Constant Flux Constant Sedimentation: assumes both constant flux
#           and constant sedimentation rate. Fits a log-linear regression.
#
# HOW serac READS DATA
# --------------------
# serac reads from a tab-delimited .txt file stored at:
#   <working_dir>/Cores/<core_name>/<core_name>.txt
#
# Required columns (exact names, matching serac_example_ALO09P12):
#   depth_top     — top of slice (mm by default; we set input_depth_mm=FALSE for cm)
#   depth_bottom  — bottom of slice (mm or cm, same as above)
#   density       — dry bulk density (g/cm³)
#   Pbex          — EXCESS (unsupported) Pb-210 activity (Bq/kg)
#   Pbex_er       — 1 SD counting error on Pbex (Bq/kg)
#   Cs            — Cs-137 activity (Bq/kg); NA if not measured
#   Cs_er         — 1 SD error on Cs-137; NA if not measured
#   Am            — Am-241 activity (Bq/kg); NA if not measured
#   Am_er         — 1 SD error on Am-241; NA if not measured
#
# UNIT CONVERSION
# ---------------
# serac expects activity in Bq/kg. Lab data is often in DPM/g.
# Conversion: 1 DPM/g = 16.667 Bq/kg
#
# STEPS
# -----
# 1. Install / load packages
# 2. Load and prepare data
# 3. Estimate supported (background) Pb-210
# 4. Build serac input file and run models
# 5. Plot and compare outputs
# =============================================================================


# =============================================================================
# STEP 1 — Install and load packages
# =============================================================================

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("serac",    quietly = TRUE)) devtools::install_github("rosalieb/serac")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr",    quietly = TRUE)) install.packages("tidyr")

library(serac)
library(dplyr)
library(ggplot2)
library(tidyr)

DPM_G_TO_BQ_KG <- 16.667   # unit conversion factor


# =============================================================================
# STEP 2 — Load and prepare data
# =============================================================================

data_file <- "RERCA/data/example_pb210_data.csv"   # <-- CHANGE THIS: paste full path to your CSV, e.g.:
# data_file <- "/Users/cathaldoherty/Desktop/CoastalBC_WorkflowV1/example_pb210_data.csv"

core_label  <- "MyCore"   # <-- short label, no spaces — used for folder and file names
coring_year <- 2024       # <-- calendar year core was collected (AD)

core_raw <- read.csv(data_file, na.strings = c("NA", "")) %>%
  mutate(
    depth_cm           = (upper_depth_cm + lower_depth_cm) / 2,
    slice_thickness_cm = lower_depth_cm - upper_depth_cm
  )

# Interpolate missing bulk density
if (any(is.na(core_raw$dry_bulk_density_g_cm3))) {
  cat(sprintf("Interpolating %d missing bulk density value(s).\n",
              sum(is.na(core_raw$dry_bulk_density_g_cm3))))
  core_raw <- core_raw %>%
    mutate(dry_bulk_density_g_cm3 = approx(
      x    = depth_cm[!is.na(dry_bulk_density_g_cm3)],
      y    = dry_bulk_density_g_cm3[!is.na(dry_bulk_density_g_cm3)],
      xout = depth_cm, rule = 2)$y)
}

cat("Data loaded. Rows:", nrow(core_raw), "\n")
print(head(core_raw %>% select(upper_depth_cm, lower_depth_cm,
                                dry_bulk_density_g_cm3, pb210_total_dpm_g)))


# =============================================================================
# STEP 3 — Estimate supported (background) Pb-210
# =============================================================================

has_ra226 <- "ra226_supported_dpm_g" %in% names(core_raw) &&
             sum(!is.na(core_raw$ra226_supported_dpm_g)) >= 3

if (has_ra226) {
  core_raw <- core_raw %>%
    mutate(pb210_supported     = ra226_supported_dpm_g,
           pb210_supported_err = ifelse("ra226_error_dpm_g" %in% names(.),
                                        ra226_error_dpm_g, ra226_supported_dpm_g * 0.05))
  cat("Supported Pb-210: from Ra-226 measurements.\n")
} else {
  cat("No Ra-226 data. Estimating supported activity from deep-core tail.\n\n")
  cat("Choose method:\n")
  cat("  1 = Automatic (mean of deepest 20% of measured sections)\n")
  cat("  2 = Manual (enter your own value)\n")
  method_choice <- readline("Enter 1 or 2: ")

  if (trimws(method_choice) == "2") {
    bg_value <- as.numeric(readline("Supported Pb-210 activity (DPM/g): "))
    bg_error <- as.numeric(readline("Uncertainty (1 SD, DPM/g): "))
  } else {
    measured   <- core_raw %>% filter(!is.na(pb210_total_dpm_g))
    n_tail     <- max(2, round(nrow(measured) * 0.2))
    tail_data  <- tail(measured, n_tail)
    bg_value   <- mean(tail_data$pb210_total_dpm_g)
    bg_error   <- sd(tail_data$pb210_total_dpm_g)
    if (is.na(bg_error) || bg_error == 0) bg_error <- bg_value * 0.1
    cat(sprintf("Automatic background: %.3f ± %.3f DPM/g (%d sections)\n",
                bg_value, bg_error, n_tail))
  }

  core_raw <- core_raw %>%
    mutate(pb210_supported     = bg_value,
           pb210_supported_err = bg_error)
}

# Calculate unsupported (excess) Pb-210 and convert to Bq/kg for serac
core_raw <- core_raw %>%
  mutate(
    pb210_unsupported_dpm_g  = pb210_total_dpm_g - pb210_supported,
    pb210_unsupported_err_dpm = sqrt(pb210_error_dpm_g^2 + pb210_supported_err^2),
    Pbex_Bq_kg               = pb210_unsupported_dpm_g  * DPM_G_TO_BQ_KG,
    Pbex_er_Bq_kg            = pb210_unsupported_err_dpm * DPM_G_TO_BQ_KG
  )

cat("\n--- Unsupported Pb-210 profile (DPM/g → Bq/kg) ---\n")
print(core_raw %>%
  filter(!is.na(pb210_total_dpm_g)) %>%
  select(depth_cm, pb210_unsupported_dpm_g, Pbex_Bq_kg, Pbex_er_Bq_kg))


# =============================================================================
# STEP 4 — Build serac input file and run models
# =============================================================================
# The input file must match the format of serac_example_ALO09P12 exactly.
# serac reads ALL rows from the file — include sections with NA activity.

# Build input with all rows (serac handles NAs internally)
serac_input <- core_raw %>%
  transmute(
    depth_top    = upper_depth_cm,          # in cm; we pass input_depth_mm=FALSE
    depth_bottom = lower_depth_cm,
    density      = dry_bulk_density_g_cm3,
    Pbex         = Pbex_Bq_kg,
    Pbex_er      = Pbex_er_Bq_kg,
    Cs           = if ("cs137_dpm_g"       %in% names(core_raw)) cs137_dpm_g       * DPM_G_TO_BQ_KG else NA_real_,
    Cs_er        = if ("cs137_error_dpm_g" %in% names(core_raw)) cs137_error_dpm_g * DPM_G_TO_BQ_KG else NA_real_,
    Am           = NA_real_,
    Am_er        = NA_real_
  )

# Create folder structure and write the input file
serac_dir <- file.path("Cores", core_label)
dir.create(serac_dir, showWarnings = FALSE, recursive = TRUE)
serac_txt <- file.path(serac_dir, paste0(core_label, ".txt"))
write.table(serac_input, serac_txt, row.names = FALSE, sep = "\t", na = "")

cat(sprintf("\nserac input written to: %s\n", serac_txt))
cat("Columns:", paste(names(serac_input), collapse = ", "), "\n")
cat("First rows:\n"); print(head(serac_input))

cat("\n--- Running serac (CRS + CFCS) ---\n")
cat("serac opens plots interactively — close each window to continue.\n\n")

serac_result <- tryCatch(
  serac(
    name           = core_label,
    coring_yr      = coring_year,
    model          = c("CRS", "CFCS"),
    input_depth_mm = FALSE,      # depths in our file are in cm, not mm
    dmax           = max(core_raw$lower_depth_cm, na.rm = TRUE),
    plotpdf        = FALSE,
    plottiff       = FALSE,
    preview        = FALSE,
    save_code      = FALSE
  ),
  error = function(e) {
    cat("serac() error:", conditionMessage(e), "\n")
    cat("Input file:    ", serac_txt, "\n")
    cat("Working dir:   ", getwd(), "\n")
    NULL
  }
)

cat("\nserac run complete.\n")


# =============================================================================
# STEP 5 — Extract and plot age-depth results
# =============================================================================

if (!is.null(serac_result)) {
  cat("Output names:", paste(names(serac_result), collapse = ", "), "\n")
  str(serac_result, max.level = 2)

  # Extract age tables — serac returns results in a named list
  extract_ages <- function(res, model_name) {
    if (is.null(res)) return(NULL)
    age_col   <- intersect(c("age", "Age", "mean_age", "age_mean"),     names(res))[1]
    lo_col    <- intersect(c("age_min", "age.min", "lwr", "age_min_95"), names(res))[1]
    hi_col    <- intersect(c("age_max", "age.max", "upr", "age_max_95"), names(res))[1]
    depth_col <- intersect(c("depth", "Depth", "depth_top",
                              "depth_top_cm", "depth_mid"),              names(res))[1]
    if (is.na(age_col) || is.na(depth_col)) {
      cat("Column names in", model_name, "output:", paste(names(res), collapse = ", "), "\n")
      return(NULL)
    }
    data.frame(
      depth  = res[[depth_col]],
      age    = res[[age_col]],
      age_lo = if (!is.na(lo_col)) res[[lo_col]] else NA_real_,
      age_hi = if (!is.na(hi_col)) res[[hi_col]] else NA_real_,
      model  = model_name
    )
  }

  # serac returns ages in AD calendar years and depths in mm.
  # We extract directly using the confirmed column names from the str() output.
  cfcs_raw <- serac_result[["CFCS age-depth model interpolated"]]
  crs_raw  <- serac_result[["CRS age-depth model interpolated"]]

  cfcs_ages <- if (!is.null(cfcs_raw)) {
    cfcs_raw %>%
      filter(!is.na(BestAD)) %>%
      transmute(
        depth  = depth_avg_mm / 10,              # mm → cm
        age    = coring_year - BestAD,           # AD → years before coring
        age_lo = coring_year - MaxAD,            # note: MaxAD = youngest = smallest age_BP
        age_hi = coring_year - MinAD,            # MinAD = oldest = largest age_BP
        model  = "CFCS"
      )
  }

  crs_ages <- if (!is.null(crs_raw)) {
    crs_raw %>%
      filter(!is.na(BestAD_CRS)) %>%
      transmute(
        depth  = depth_top_mm / 10,
        age    = coring_year - BestAD_CRS,
        age_lo = coring_year - MaxAD_CRS,
        age_hi = coring_year - MinAD_CRS,
        model  = "CRS"
      )
  }

  ages_all <- bind_rows(cfcs_ages, crs_ages)

  if (!is.null(ages_all) && nrow(ages_all) > 0) {
    cat("\n=== serac Age-Depth Results (years before coring) ===\n")
    print(ages_all %>% mutate(across(where(is.numeric), \(x) round(x, 1))))

    model_colours <- c(CRS = "#2c7bb6", CFCS = "#1a9641")

    p_compare <- ggplot(ages_all, aes(x = age, y = depth, colour = model)) +
      geom_ribbon(aes(xmin = age_lo, xmax = age_hi, fill = model),
                  alpha = 0.2, colour = NA) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_y_reverse(name = "Depth (cm)") +
      scale_x_continuous(name = "Age (years before coring)") +
      scale_colour_manual(values = model_colours) +
      scale_fill_manual(values = model_colours, guide = "none") +
      labs(title    = "serac: CRS vs CFCS Age-Depth Models",
           subtitle = "Shaded bands = 95% uncertainty",
           colour   = "Model") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    print(p_compare)

    # Also plot in AD calendar years
    p_ad <- ggplot(ages_all, aes(x = coring_year - age, y = depth, colour = model)) +
      geom_ribbon(aes(xmin = coring_year - age_hi, xmax = coring_year - age_lo,
                      fill = model), alpha = 0.2, colour = NA) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_y_reverse(name = "Depth (cm)") +
      scale_x_continuous(name = "Calendar year (AD)") +
      scale_colour_manual(values = model_colours) +
      scale_fill_manual(values = model_colours, guide = "none") +
      labs(title  = "serac: CRS vs CFCS — Calendar Year Scale",
           colour = "Model") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    print(p_ad)

    # Print CFCS sedimentation rate
    cfcs_sar <- serac_result[["CFCS sediment accumulation rate"]]
    if (!is.null(cfcs_sar)) {
      cat(sprintf("\nCFCS sedimentation rate: %.3f mm/yr (+/- %.3f mm/yr)  R² = %.3f\n",
                  abs(cfcs_sar$SAR_mm.yr.1), abs(cfcs_sar$error_mm.yr.1), cfcs_sar$R2))
    }

    for (mod in c("CRS", "CFCS")) {
      sub <- ages_all %>% filter(model == mod)
      if (nrow(sub) > 0)
        cat(sprintf("%s: datable range 0–%.1f cm  |  %d–%d AD  |  0–%.0f yr before coring\n",
                    mod,
                    max(sub$depth, na.rm = TRUE),
                    round(min(coring_year - sub$age, na.rm = TRUE)),
                    coring_year,
                    max(sub$age, na.rm = TRUE)))
    }
  } else {
    cat("No ages extracted — check str() output above.\n")
  }
}

# --- Activity profile for reference ---
act_data <- core_raw %>% filter(!is.na(pb210_total_dpm_g))

p_activity <- ggplot(act_data, aes(y = depth_cm)) +
  geom_point(aes(x = pb210_total_dpm_g * DPM_G_TO_BQ_KG), colour = "#2c7bb6", size = 2.5) +
  geom_errorbarh(aes(xmin = (pb210_total_dpm_g - pb210_error_dpm_g) * DPM_G_TO_BQ_KG,
                     xmax = (pb210_total_dpm_g + pb210_error_dpm_g) * DPM_G_TO_BQ_KG),
                 colour = "#2c7bb6", height = 0.3) +
  geom_point(aes(x = Pbex_Bq_kg), colour = "#d7191c", size = 2.5, shape = 17) +
  geom_vline(xintercept = bg_value * DPM_G_TO_BQ_KG,
             linetype = "dashed", colour = "grey50") +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Activity (Bq/kg)") +
  labs(title    = "Pb-210 Activity Profile",
       subtitle = "Blue = total; red = unsupported; dashed = background") +
  theme_bw(base_size = 12)
print(p_activity)

# =============================================================================
# SAVE OUTPUTS for script 05_summary.R
# =============================================================================
dir.create("RERCA/output", showWarnings = FALSE, recursive = TRUE)

if (exists("ages_all") && !is.null(ages_all) && nrow(ages_all) > 0) {
  serac_out <- ages_all %>%
    rename(depth_cm = depth, age_min_yr = age_lo, age_max_yr = age_hi) %>%
    mutate(model = paste0("serac_", model))
  write.csv(serac_out, "RERCA/output/03_serac_ages.csv", row.names = FALSE)

  if (exists("cfcs_sar") && !is.null(cfcs_sar)) {
    stats_out <- data.frame(
      model            = c("serac_CRS", "serac_CFCS"),
      method           = c("CRS (serac)", "CFCS (serac)"),
      package          = "serac",
      cfcs_sar_mm_yr   = c(NA, abs(cfcs_sar$SAR_mm.yr.1)),
      cfcs_sar_err     = c(NA, abs(cfcs_sar$error_mm.yr.1)),
      cfcs_r2          = c(NA, cfcs_sar$R2),
      stringsAsFactors = FALSE
    )
    write.csv(stats_out, "RERCA/output/03_serac_stats.csv", row.names = FALSE)
  }
  cat("Outputs saved to RERCA/output/\n")
}

cat("\nScript 03 complete.\n")
cat("Next: run 04_bayesian_pb210.R for the manual Bayesian approach.\n")
