library(glmmTMB)
library(dplyr)
library(ggplot2)
library(future.apply)

# ==============================================================================
# 1. PARAMETER EXTRACTION & HIGH-PRECISION CALIBRATION
# ==============================================================================
# Load data (Adjust path as needed)
# trout <- read.csv("G:/My Drive/IBBEAM/Power analysis/Raw data/Sportfish/trout.csv") 
trout <- read.csv("C:/Users/nal102/OneDrive - University of Miami/Desktop/Habitat ecology/Sportfish project/CSV/trout_counts_&_env.csv") %>%
  filter(year <= 2017, month %in% 5:11) %>%
  mutate(sg_cover = mean_pc_Thalassia + mean_pc_Halodule + mean_pc_Syringodium,
         zone = as.factor(zone),
         month = as.factor(month))



# Apply Sum Coding for unbiased baseline estimation
contrasts(trout$zone) <- contr.sum(levels(trout$zone))

mod_pa_pilot <- glmmTMB(occurrence ~ zone + 
                          scale(temp, scale = FALSE) +
                          scale(sal, scale = FALSE) + 
                          splines::ns(scale(sg_cover, scale = FALSE), 3) +
                          (1|station) + (1|year) + (1|month), 
                        data = trout, family = binomial, REML = TRUE)

# Extract Parameters
fe        <- fixef(mod_pa_pilot)$cond
zone_effs <- fe[grep("zone", names(fe))]
re_var    <- VarCorr(mod_pa_pilot)$cond
site_sd   <- attr(re_var$station, "stddev")
year_sd   <- attr(re_var$year, "stddev")
month_sd  <- attr(re_var$month, "stddev")
obs_mean  <- mean(trout$occurrence, na.rm = TRUE)
total_sd  <- sqrt(site_sd^2 + year_sd^2 + month_sd^2)

# Reconstitute Zone Offsets
levs <- levels(trout$zone)
z_map <- setNames(rep(0, length(levs)), levs)
for(z_name in names(zone_effs)) {
  clean_name <- gsub("zone", "", z_name)
  z_map[clean_name] <- zone_effs[z_name]
}
z_map[levs[length(levs)]] <- -sum(zone_effs)

# --- CALIBRATION ---
find_intercept <- function(candidate_logit, z_lookup, target, sd_val, levels_vec) {
  n_test <- 100000
  test_df <- data.frame(zone = factor(rep(levels_vec, length.out = n_test)))
  offsets <- as.numeric(z_lookup[as.character(test_df$zone)])
  noise   <- rnorm(n_test, 0, sd_val)
  sim_p <- mean(plogis(candidate_logit + offsets + noise))
  return(sim_p - target)
}

calib_fit <- uniroot(find_intercept, c(-20, 2), 
                     z_lookup = z_map, target = obs_mean, 
                     sd_val = total_sd, levels_vec = levs)
base_logit <- calib_fit$root

# Baseline intercept starts here:
cat("Target Presence:     ", round(obs_mean * 100, 2), "%\n")
cat("Calibrated Intercept:", round(base_logit, 3), "\n")
cat("==============================================\n\n")

# ==============================================================================
# 2. CORE SIMULATION FUNCTIONS (OPTIMIZED)
# ==============================================================================

plan(multisession, workers = parallel::detectCores() - 1)

sim_trout_fast <- function(skeleton, n_sites, n_years, site_sd, year_sd, month_sd) {
  
  # A. GENERATIVE PHASE (Vectorized Random Effects)
  # Uses indexing instead of left_join for massive speed boost
  u_site  <- rnorm(n_sites, 0, site_sd)[as.numeric(skeleton$site)]
  u_year  <- rnorm(n_years, 0, year_sd)[skeleton$year + 1]
  
  # Monthly random effects (unique per year-month combination)
  u_month_vec <- rnorm(n_years * 7, 0, month_sd)
  u_month     <- u_month_vec[as.numeric(interaction(skeleton$month_id, skeleton$year))]
  
  # Combine with pre-calculated fixed effects
  logit_p <- skeleton$lp_fixed + u_site + u_year + u_month
  skeleton$presence <- rbinom(nrow(skeleton), 1, plogis(logit_p))
  
  # B. AGGREGATION PHASE
  df_agg <- skeleton %>% 
    group_by(site, year, zone) %>%
    summarise(positives = sum(presence), 
              total = n(), 
              negatives = total - positives, 
              .groups = "drop")
  
  # C. ESTIMATION PHASE
  # Using betabinomial to capture year/month noise without a random year effect
  mod <- try(glmmTMB(cbind(positives, negatives) ~ year + zone + (1|site), 
                     data = df_agg, 
                     family = betabinomial,
                     control = glmmTMBControl(optimizer = nlminb, 
                                              optCtrl = list(iter.max = 400, eval.max = 400),
                                              eigval_check = FALSE)), silent = TRUE)
  
  if(!inherits(mod, "try-error")) {
    summ <- try(summary(mod)$coefficients$cond, silent = TRUE)
    if(!inherits(summ, "try-error") && "year" %in% rownames(summ)) {
      p_val <- summ["year", "Pr(>|z|)"]
      est   <- summ["year", "Estimate"]
      # Two-tailed significance check at alpha = 0.05
      if(!is.na(p_val) && !is.na(est)) return((p_val / 2) < 0.05 && est > 0)
    }
  }
  return(FALSE)
}

run_parallel_power <- function(n_sites, total_change, n_years, n_sims, ... ) {
  
  # Pre-calculate Trend Logic
  p_start <- obs_mean 
  p_end   <- p_start * (total_change)
  if(p_end >= 0.99) p_end <- 0.99
  
  total_logit_change <- qlogis(p_end) - qlogis(p_start)
  ann_slope <- total_logit_change / (n_years - 1)
  
  # Create Skeleton ONCE per site-count iteration
  skeleton <- expand.grid(site = factor(1:n_sites), 
                          year = 0:(n_years - 1), 
                          month_id = factor(5:11))
  
  site_map <- data.frame(site = factor(1:n_sites), 
                         zone = factor(rep(levels(trout$zone), length.out = n_sites)))
  skeleton <- left_join(skeleton, site_map, by = "site")
  
  # PRE-CALCULATE FIXED LINEAR PREDICTOR
  skeleton$lp_fixed <- base_logit + 
    as.numeric(z_map[as.character(skeleton$zone)]) + 
    (ann_slope * skeleton$year)
  
  cat("Running for n_sites =", n_sites, "...\n")
  
  results <- future_lapply(1:n_sims, function(x) {
    sim_trout_fast(skeleton, n_sites, n_years, ...)
  }, future.seed = TRUE)
  
  return(mean(as.logical(unlist(results))) * 100)
}

# ==============================================================================
# 3. EXECUTION
# ==============================================================================

site_counts <- seq(80, 160, by = 20)
# site_counts <- 112

plot_data_pa <- data.frame(
  n_sites = site_counts,
  power = sapply(site_counts, function(n) {
    res <- run_parallel_power(
      n_sims = 1000, 
      n_sites = n, 
      total_change = 1.3, 
      n_years = 5, 
      site_sd = site_sd, 
      year_sd = year_sd, 
      month_sd = month_sd
    )
    return(as.numeric(res))
  })
)

length(unique(trout$station))
table(trout$month, trout$year) # 81

total_changes <- seq(2.0, 4.0, by=0.2)  # your list

plot_factor_inc <- data.frame(
  total_change = total_changes,
  power = sapply(total_changes, function(tc) {
    res <- run_parallel_power(
      n_sims = 1000, 
      n_sites = 81,   # fixed value
      total_change = tc, 
      n_years = 5, 
      site_sd = site_sd, 
      year_sd = year_sd, 
      month_sd = month_sd
    )
    as.numeric(res)
  })
)

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

# Relative change (compared to starting point)
ggplot(plot_data_pa, aes(x = n_sites, y = power)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "black", size = 4) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(title = "Seatrout Presence-Absence Power Analysis (Beta-Binomial)",
       subtitle = "Detecting +30% Δ from 6.4 prob. of occurrence in 5 years",
       x = "Number of Sites (per month)", y = "Power (%)") +
  theme_minimal(base_size = 14)

# expected occurrences by zone
emmeans::emmeans(mod_pa_pilot, "zone", type = "response",
                 bias.adj = TRUE, sigma = sqrt(sum(unlist(VarCorr(mod_pa_pilot)))))

# (Factor - 1) * 100 = relative % change
# A 6.4% percent change of seeing a trout at the average site
ggplot(plot_data_pa, aes(x = total_change, y = power)) +
  geom_line(color = "#16a085", linewidth = 1.2) +
  geom_point(color = "#2c3e50", size = 3) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(
    title = "Seatrout Presence-Absence Power Analysis (Beta-Binomial)",
    subtitle = "Detecting Δ from 6.4 prob. of occurrence in 5 years",
    x = "Multi. factor of baseline occurrence (at 81 sites per month)", 
    y = "Power (%)"
  ) +
  theme_minimal(base_size = 14)

emmeans::emmeans(mod_pa_pilot, ~ 1, 
                 type = "response",
                 bias.adj = TRUE, # Jensen's Adjustment
                 sigma = sqrt(sum(unlist(VarCorr(mod_pa_pilot)))))

# 1        prob     SE  df asymp.LCL asymp.UCL
# overall 0.064 0.0183 Inf     0.036      0.11
# 
# Results are averaged over the levels of: zone 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 
# Bias adjustment applied based on sigma = 1.1876 

#' “Across all sites, months, years, and zones, accounting for random variability, 
#' `the expected probability of trout occurrence is ~6.4%.”`

plot_data_pa
#    total_change power
# 1           2.0  58.5
# 2           2.2  58.9
# 3           2.4  65.8
# 4           2.6  65.6
# 5           2.8  74.3
# 6           3.0  72.7
# 7           3.2  77.0
# 8           3.4  79.5
# 9           3.6  79.4
# 10          3.8  82.5
# 11          4.0  83.3


