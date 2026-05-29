# =============================================================================
# RERCA | Script 04: Manual Bayesian Pb-210 Chronology
# Reference: Dunnington (2019) — https://dewey.dunnington.ca/post/2019/
#            doing-bayesian-lead-210-interpretation/
# =============================================================================
#
# OVERVIEW
# --------
# This script implements the manual Bayesian approach described by Dunnington
# (2019). Rather than a black-box package, the model is built explicitly in
# Stan via the {rstan} or {cmdstanr} interface, giving you full control over
# priors, the likelihood, and how uncertainty propagates.
#
# The model simultaneously estimates:
#   - The supported (background) Pb-210 activity (phi)
#   - The initial surface activity (A0)
#   - The mean sedimentation rate and its variability
#   - A posterior age for every dated depth
#
# STEPS
# -----
# 1. Install / load packages
# 2. Load and prepare data
# 3. Define and compile the Stan model
# 4. Run MCMC sampling
# 5. Inspect diagnostics (Rhat, n_eff, trace plots)
# 6. Extract and plot posterior age-depth estimates
# 7. Compare with CRS (optional)
# =============================================================================


# =============================================================================
# STEP 1 — Install and load packages
# =============================================================================
# rstan requires a C++ toolchain. If you have not set this up before, visit:
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

if (!requireNamespace("rstan",   quietly = TRUE)) install.packages("rstan")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("bayesplot", quietly = TRUE)) install.packages("bayesplot")

library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bayesplot)

# Stan-wide options — use all available CPU cores for parallel chains
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# =============================================================================
# STEP 2 — Load and prepare data
# =============================================================================

data_file <- "RERCA/data/example_pb210_data.csv"   # <-- CHANGE THIS: paste the full path to your CSV here, e.g.:

core_raw <- read.csv(data_file, na.strings = c("NA", "")) %>%
  mutate(
    depth_cm           = (upper_depth_cm + lower_depth_cm) / 2,
    slice_thickness_cm = lower_depth_cm - upper_depth_cm
  )

# Interpolate missing bulk density
if (any(is.na(core_raw$dry_bulk_density_g_cm3))) {
  core_raw <- core_raw %>%
    mutate(dry_bulk_density_g_cm3 = approx(
      x    = depth_cm[!is.na(dry_bulk_density_g_cm3)],
      y    = dry_bulk_density_g_cm3[!is.na(dry_bulk_density_g_cm3)],
      xout = depth_cm, rule = 2)$y)
}

core_raw <- core_raw %>%
  mutate(
    mass_depth_g_cm2 = cumsum(dry_bulk_density_g_cm3 * slice_thickness_cm)
  )

# Keep only sections with measured Pb-210
measured <- core_raw %>% filter(!is.na(pb210_total_dpm_g))

cat(sprintf("Sections with Pb-210 measurements: %d\n", nrow(measured)))
cat("Mass depth range: 0 –", round(max(measured$mass_depth_g_cm2), 2), "g/cm²\n")

# Inspect the activity profile before modelling
p_raw <- ggplot(measured, aes(x = pb210_total_dpm_g, y = depth_cm)) +
  geom_point(size = 2.5, colour = "#2c7bb6") +
  geom_errorbarh(aes(xmin = pb210_total_dpm_g - pb210_error_dpm_g,
                     xmax = pb210_total_dpm_g + pb210_error_dpm_g),
                 height = 0.3, colour = "#2c7bb6") +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Total Pb-210 (DPM / g dry wt)") +
  labs(title = "Raw Pb-210 Activity Profile",
       subtitle = "Review before setting priors") +
  theme_bw(base_size = 12)
print(p_raw)


# =============================================================================
# STEP 3 — Define the Stan model
# =============================================================================
# The model assumes:
#   A(x) = phi + A0 * exp(-lambda * t(x))
#   where:
#     A(x)   = total activity at cumulative mass depth x
#     phi    = supported (background) activity
#     A0     = initial unsupported surface activity
#     lambda = Pb-210 decay constant (0.03114 yr-1, fixed)
#     t(x)   = age at depth x = x / (rho * omega), integrated with variable rate
#
# For variable sedimentation, age at depth i is modelled as the sum of
# time increments delta_t[i] ~ lognormal(log(mean_acc), sigma_acc)

lambda_pb210 <- 0.03114  # yr-1 (half-life = 22.26 yr)

stan_model_code <- "
data {
  int<lower=1> N;                    // number of measured sections
  vector[N] activity;                // total Pb-210 activity (DPM/g)
  vector<lower=0>[N] activity_err;   // 1 SD measurement error
  vector<lower=0>[N] mass_depth;     // cumulative mass depth (g/cm2)
  real<lower=0> lambda;              // decay constant (yr-1)
}

parameters {
  real<lower=0> phi;          // supported activity (DPM/g)
  real<lower=0> A0;           // initial unsupported surface activity (DPM/g)
  real<lower=0> mean_acc;     // mean accumulation rate (yr / (g/cm2))
  real<lower=0> sigma_acc;    // variability in accumulation rate (log scale)
  vector<lower=0>[N] delta_t; // time increment for each section (yr)
}

transformed parameters {
  vector[N] age;             // age at each measured depth (yr)
  vector[N] pred_activity;   // predicted total activity
  age[1] = delta_t[1];
  for (i in 2:N)
    age[i] = age[i-1] + delta_t[i];
  for (i in 1:N)
    pred_activity[i] = phi + A0 * exp(-lambda * age[i]);
}

model {
  // Priors — adjust these based on your knowledge of the system
  phi      ~ normal(0.5, 0.5);         // supported activity: ~0.5 DPM/g
  A0       ~ normal(2.0, 2.0);         // initial surface activity
  mean_acc ~ normal(5.0, 5.0);         // mean time per unit mass depth (yr / g/cm2)
  sigma_acc ~ exponential(1.0);

  // Each time increment drawn from a lognormal centred on mean_acc
  for (i in 1:N)
    delta_t[i] ~ lognormal(log(mean_acc), sigma_acc);

  // Likelihood: observed activity is Gaussian around predicted
  for (i in 1:N)
    activity[i] ~ normal(pred_activity[i], activity_err[i]);
}
"

cat("Compiling Stan model...\n")
stan_model_obj <- stan_model(model_code = stan_model_code, model_name = "pb210_bayesian")
cat("Compilation complete.\n")


# =============================================================================
# STEP 4 — Run MCMC sampling
# =============================================================================
# PRIOR GUIDANCE:
#   phi      — supported activity. If you can see the deep-core tail plateauing
#              in the plot above, use that value ± uncertainty as your prior mean/sd.
#   A0       — surface activity. Read off the approximate value at 0 cm from the
#              plot. Default prior is broad (Normal(2, 2)) — adjust if needed.
#   mean_acc — prior mean for yr per g/cm². For a typical lake, 5–20 yr/(g/cm²)
#              is a reasonable starting range.
#
# You can tighten priors by editing the `model {}` block above.

stan_data <- list(
  N            = nrow(measured),
  activity     = measured$pb210_total_dpm_g,
  activity_err = measured$pb210_error_dpm_g,
  mass_depth   = measured$mass_depth_g_cm2,
  lambda       = lambda_pb210
)

cat("\nRunning MCMC — 4 chains × 2000 iterations (1000 warmup)...\n")
fit <- sampling(
  stan_model_obj,
  data    = stan_data,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\nSampling complete.\n")


# =============================================================================
# STEP 5 — Diagnostics
# =============================================================================

cat("\n=== Model Diagnostics ===\n")
print(summary(fit, pars = c("phi", "A0", "mean_acc", "sigma_acc"))$summary)

# Rhat should be close to 1.00 (< 1.01 is excellent, < 1.05 is acceptable)
# n_eff should be > 400 per parameter

# Check for divergences
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat(sprintf("\nDivergent transitions: %d\n", n_divergent))
if (n_divergent > 0) {
  cat("Consider increasing adapt_delta (currently 0.95) or adjusting priors.\n")
} else {
  cat("No divergences — sampler converged well.\n")
}

# Trace plots for key parameters
posterior_array <- as.array(fit)
p_trace <- mcmc_trace(posterior_array, pars = c("phi", "A0", "mean_acc", "sigma_acc"))
print(p_trace)

# Posterior distributions
p_dens <- mcmc_dens_overlay(posterior_array, pars = c("phi", "A0", "mean_acc"))
print(p_dens)


# =============================================================================
# STEP 6 — Extract and plot age-depth estimates
# =============================================================================

age_samples <- rstan::extract(fit, pars = "age")$age
# age_samples is a matrix: rows = MCMC samples, columns = depth sections

age_summary <- data.frame(
  depth_cm   = measured$depth_cm,
  age_mean   = apply(age_samples, 2, mean),
  age_median = apply(age_samples, 2, median),
  age_lo_95  = apply(age_samples, 2, quantile, probs = 0.025),
  age_hi_95  = apply(age_samples, 2, quantile, probs = 0.975),
  age_lo_50  = apply(age_samples, 2, quantile, probs = 0.25),
  age_hi_50  = apply(age_samples, 2, quantile, probs = 0.75)
)

cat("\n=== Bayesian Posterior Age Estimates ===\n")
print(age_summary %>% mutate(across(where(is.numeric), \(x) round(x, 1))))

# Age-depth plot
p_ages <- ggplot(age_summary, aes(x = age_mean, y = depth_cm)) +
  geom_ribbon(aes(xmin = age_lo_95, xmax = age_hi_95),
              fill = "#d7191c", alpha = 0.2) +
  geom_ribbon(aes(xmin = age_lo_50, xmax = age_hi_50),
              fill = "#d7191c", alpha = 0.4) +
  geom_line(colour = "#d7191c", linewidth = 1) +
  geom_point(colour = "#d7191c", size = 2.5) +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Age (years before coring)") +
  labs(title   = "Bayesian Pb-210 Age-Depth Model (Dunnington 2019)",
       subtitle = "Dark band = 50% CI; light band = 95% CI") +
  theme_bw(base_size = 12)
print(p_ages)

# Posterior predictive check — do predicted activities match observations?
phi_samples   <- rstan::extract(fit, pars = "phi")$phi
A0_samples    <- rstan::extract(fit, pars = "A0")$A0

pred_check <- data.frame(
  depth_cm      = measured$depth_cm,
  obs_activity  = measured$pb210_total_dpm_g,
  obs_err       = measured$pb210_error_dpm_g,
  pred_mean     = apply(rstan::extract(fit, pars = "pred_activity")$pred_activity, 2, mean),
  pred_lo       = apply(rstan::extract(fit, pars = "pred_activity")$pred_activity, 2, quantile, 0.025),
  pred_hi       = apply(rstan::extract(fit, pars = "pred_activity")$pred_activity, 2, quantile, 0.975)
)

p_ppc <- ggplot(pred_check, aes(y = depth_cm)) +
  geom_ribbon(aes(xmin = pred_lo, xmax = pred_hi), fill = "#fdae61", alpha = 0.4) +
  geom_line(aes(x = pred_mean), colour = "#fdae61", linewidth = 1) +
  geom_point(aes(x = obs_activity), colour = "#2c7bb6", size = 2.5) +
  geom_errorbarh(aes(xmin = obs_activity - obs_err, xmax = obs_activity + obs_err),
                 colour = "#2c7bb6", height = 0.3) +
  scale_y_reverse(name = "Depth (cm)") +
  scale_x_continuous(name = "Pb-210 Activity (DPM / g dry wt)") +
  labs(title   = "Posterior Predictive Check",
       subtitle = "Blue = observed ± 1 SD; orange = model 95% credible interval") +
  theme_bw(base_size = 12)
print(p_ppc)

cat(sprintf(
  "\nPosterior supported activity (phi): mean = %.3f, 95%% CI = %.3f – %.3f DPM/g\n",
  mean(phi_samples),
  quantile(phi_samples, 0.025),
  quantile(phi_samples, 0.975)
))
cat(sprintf(
  "Posterior surface activity  (A0):  mean = %.3f, 95%% CI = %.3f – %.3f DPM/g\n",
  mean(A0_samples),
  quantile(A0_samples, 0.025),
  quantile(A0_samples, 0.975)
))


# =============================================================================
# STEP 7 — Compare all methods (optional)
# =============================================================================
# If you have run scripts 01, 02, and 03 and saved their outputs, load them
# here to make a direct comparison plot.

cat("\nScript 04 complete.\n")
cat("You have now run all four Pb-210 dating approaches.\n")
cat("Compare the age-depth plots to assess consistency between methods.\n")
cat("Consistent results across CRS, rplum, serac, and Bayesian models increase confidence.\n")
