
#' power_for_n_after
#'
#' Estimate power for a paired t-test on the average of a given number of
#' measurements at a fixed number of sites before and after a changepoint.
#'
#' @param S Number of sites
#' @param nB  Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param delta Hypothesized mean change
#' @param sd_w Within-site standard deviation
#' @param sd_d Between-site standard deviation of true changes
#' @param alpha Significance level, defaults to 0.05
#' @param nsim Number of simulations to run
#' @param seed Random seed for reproducibility
#'
#' @returns Estimated power for the given number of after measurements
#' @export
#'
power_for_n_after <- function(S, nB, nA,
                              delta, sd_w, sd_d = 0,
                              alpha = 0.05, nsim = 2000,
                              seed = 1) {
  set.seed(seed)
  pvals <- replicate(nsim, {
    true_change <- rnorm(S, mean = delta, sd = sd_d)
    yB <- matrix(rnorm(S * nB, mean = 0, sd = sd_w), nrow = S)
    yA <- matrix(rnorm(S * nA, mean = true_change, sd = sd_w), nrow = S)
    D <- rowMeans(yA) - rowMeans(yB)
    t.test(D, mu = 0)$p.value
  })
  mean(pvals < alpha)
}

#' find_n_after
#'
#' Function to find the minimum number of after measurements needed to achieve target
#' power, for a paired t-test on the before and after means.
#'
#' @param S Number of sites
#' @param nB Number of before measurements per site
#' @param delta Hypothesized mean change
#' @param sd_w Within-site standard deviation
#' @param sd_d Between-site standard deviation of true changes
#' @param target_power Desired power level to achieve
#' @param alpha Significance level, defaults to 0.05
#' @param n_grid Grid of after measurements to evaluate
#' @param nsim Number of simulations to run
#' @param seed Random seed for reproducibility
#'
#' @returns A list containing the minimum number of after measurements needed to achieve the target power and a data frame of the power curve
#' @export
#'
find_n_after <- function(S, nB,
                         delta, sd_w, sd_d = 0,
                         target_power = 0.8,
                         alpha = 0.05,
                         n_grid = 1:50,
                         nsim = 2000, seed = 1) {
  pow <- sapply(n_grid, function(nA) {
    power_for_n_after(S = S, nB = nB, nA = nA,
                      delta = delta, sd_w = sd_w, sd_d = sd_d,
                      alpha = alpha, nsim = nsim, seed = seed)
  })
  out <- data.frame(n_after = n_grid, power = pow)
  n_star <- out$n_after[which(out$power >= target_power)[1]]
  list(n_star = n_star, curve = out)
}

#' power_for_sites
#'
#' Function to estimate power for a paired t-test on the average of a given number
#' of sites before and after a changepoint.
#'
#' @param S Number of sites
#' @param nB Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param delta Hypothesized mean change
#' @param sd_w Within-site standard deviation
#' @param sd_d Between-site standard deviation of true changes
#' @param alpha Significance level, defaults to 0.05
#' @param nsim Number of simulations to run
#' @param seed Random seed for reproducibility
#'
#' @returns Estimated power for the given number of sites
#' @export
#'
power_for_sites <- function(S, nB, nA,
                            delta, sd_w, sd_d = 0,
                            alpha = 0.05, nsim = 2000, seed = 1) {
  power_for_n_after(S = S, nB = nB, nA = nA,
                    delta = delta, sd_w = sd_w, sd_d = sd_d,
                    alpha = alpha, nsim = nsim, seed = seed)
}

#' find_min_sites
#'
#' Function to find the minimum number of sites needed to achieve target power from
#' a paired t test on the means across observations at sites before and after a changepoint.
#'
#' @param nB Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param delta Hypothesized mean change
#' @param sd_w Within-site standard deviation
#' @param sd_d Between-site standard deviation of true changes
#' @param target_power Desired power level to achieve
#' @param alpha Significance level, defaults to 0.05
#' @param S_grid Grid of site numbers to evaluate
#' @param nsim Number of simulations to run
#' @param seed Random seed for reproducibility
#'
#' @returns A list containing the minimum number of sites needed to achieve the target power and a data frame of the power curve
#' @export
#'
find_min_sites <- function(nB, nA,
                           delta, sd_w, sd_d = 0,
                           target_power = 0.8,
                           alpha = 0.05,
                           S_grid = 2:50,
                           nsim = 2000, seed = 1) {
  S_grid <- S_grid[S_grid >= 2]
  pow <- sapply(S_grid, function(S) {
    power_for_sites(S = S, nB = nB, nA = nA,
                    delta = delta, sd_w = sd_w, sd_d = sd_d,
                    alpha = alpha, nsim = nsim, seed = seed)
  })
  out <- data.frame(S = S_grid, power = pow)
  S_star <- out$S[which(out$power >= target_power)[1]]
  list(S_star = S_star, curve = out)
}

#' getSD_within
#'
#' Estimate within-site standard deviation from baseline data.
#'
#' @param baseline Data frame containing baseline measurements
#' @param siteVar Name of the site variable in baseline data
#' @param responseVar Name of the response variable in baseline data
#'
#' @returns Estimated within-site standard deviation
#' @export
#'
getSD_within <- function(baseline,siteVar="site",responseVar="y"){
  sd_within_hat <- baseline %>%
    group_by(!!sym(siteVar)) %>%
    summarise(sd_site = sd(!!sym(responseVar))) %>%
    summarise(sd_within = mean(sd_site)) %>%
    pull(sd_within)
  return(sd_within_hat)
}

#' getSD_difference
#'
#' Estimate standard deviation of difference in means before and after change for
#' paired t-test of difference in means.
#'
#' @param sd_within standard deviation within sites
#' @param nA number of after samples per site
#' @param nB number of before samples per site
#' @param sd_delta standard deviation of true changes among sites
#'
#' @returns Standard deviation of the difference in means
#'
#' @export
#'
getSD_difference <- function(sd_within, nA, nB, sd_delta){
  sd_diff <- sqrt(sd_within^2*(1/nA + 1/nB) + sd_delta^2)
  return(sd_diff)
}

#' find_min_sites_analytical
#'
#' Function to use power.t.test to calculate power for a paired t-test on the means across
#'
#' @param n_grid Number of sites grid
#' @param delta_target Target mean change
#' @param sd_diff Standard deviation of difference in means
#'
#' @returns Power for each number of sites in n_grid
#' @export
#'
find_min_sites_analytical <- function(S_grid, delta_target, sd_diff){
  S_grid <- S_grid[S_grid >=2]
  sapply(S_grid, function(n) {
    power.t.test(
      n = n,
      delta = delta_target,
      sd = sd_diff,
      sig.level = 0.05,
      type = "one.sample",
      alternative = "two.sided"
    )$power
  } )
}

#' power_for_nA
#' Function to calculate power for given number of after samples for a paire t-test
#' of the difference in means before vs. after change using power.t.test
#'
#' @param nA Number of after samples per site
#' @param S Number of sites
#' @param nB Number of before samples per site
#' @param delta Hypothesized mean change
#' @param sd_within Standard deviation within sites
#' @param sd_delta Standard deviation of true changes among sites
#' @param alpha Significance level
#'
#' @returns Power for the given number of after sample
#' @export
#'
power_for_nA <- function(nA, S, nB, delta, sd_within, sd_delta, alpha){
  sd_diff <- getSD_difference(sd_within, nA, nB, sd_delta)
  power_result <- power.t.test(n = S,
                               delta = delta,
                               sd = sd_diff,
                               sig.level = alpha,
                               power = NULL,
                               type = "paired",
                               alternative = "two.sided")
  return(power_result$power)
}

#' summarize_baseline
#'
#' Summarize baseline data to estimate within-site and between-site standard deviations, and
#' optionally log-transform the response variable to estimate these on the log scale.
#'
#' @param baseline Data frame containing baseline measurements
#' @param siteVar Name of the site variable in baseline data
#' @param responseVar Name of the response variable in baseline data
#' @param logTransform Logical indicating whether to log-transform the response variable for estimating standard deviations on the log scale
#' @param logAdd Value to add to response variable before log-transforming to avoid issues with zeros (default 0)
#'
#' @returns A data frame containing the estimated within-site and between-site standard deviations, and the number of sites and before samples per site
#' @export
#'
summarize_baseline <- function(baseline,siteVar="site",responseVar="y",
                               logTransform=FALSE,
                               logAdd=0){
  returnVal<-baseline %>%
    group_by(!!sym(siteVar)) %>%
    summarize(
      site_mean = mean(!!sym(responseVar)),
      site_sd = sd(!!sym(responseVar)),
      site_logmean = if_else(logTransform, mean(log(!!sym(responseVar)+logAdd)), NA_real_),
      site_logsd = if_else(logTransform, sd(log(!!sym(responseVar)+logAdd)), NA_real_),
      n = n()
    ) %>%
    ungroup() %>%
    summarize(sd_within = mean(site_sd),
              sd_between = sd(site_mean),
              logsd_within = if_else(logTransform, mean(site_logsd), NA_real_),
              logsd_between = if_else(logTransform, sd(site_logmean), NA_real_),
              n_sites = n(),
              nB=if_else(length(unique(n))==1,unique(n)[1],NA)
    )
  if(any(is.na(returnVal$nB))){
    warning("Number of before samples per site is not consistent across sites. nB is set to NA.")
  }
  return(returnVal)
}




#' find_min_detectable_percent
#'
#' Calculates the smallest percentage change detectable for a given sample size.
#'
#' @param S Number of sites
#' @param nB Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param sd_within Within-site standard deviation (log calculated if logTransform=TRUE)
#' @param sd_delta Between-site SD of true changes (default 0)
#' @param logTransform Logical indicating whether to calculate the detectable change on the log scale (default TRUE)
#' @param logAdd Value to add to response variable before log-transforming to avoid issues
#' @param target_power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#'
#' @returns The minimum detectable percentage change (e.g., 30 for a 30% change)
#' @export

find_min_detectable_percent <- function(S, nB, nA,
                                        sd_within,
                                        sd_delta = 0,
                                        logTransform = FALSE,
                                        logAdd=0,
                                        target_power = 0.8,
                                        alpha = 0.05) {

  # 1. Define the root-finding function for the log-delta
  power_root_func <- function(delta) {
    # Calculate SD of the difference using your provided logic
    sd_diff <- getSD_difference(sd_within, nA, nB, sd_delta)

    # Analytical power for the log-difference
    current_power <- power.t.test(
      n = S,
      delta = delta,
      sd = sd_diff,
      sig.level = alpha,
      type = "one.sample",
      alternative = "two.sided"
    )$power

    return(current_power - target_power)
  }

  # 2. Find the required log_delta (searching between 0.001 and 5)
  fit <- uniroot(power_root_func, interval = c(1e-4, 5)) #solve equation
  log_delta_required <- fit$root

  # 3. Convert log-scale change back to a percentage
  # exp(log_delta) gives the ratio (e.g., 1.3).
  # (ratio - 1) * 100 gives the percent (e.g., 30)
  if(logTransform) percent_change <- (exp(log_delta_required) - logAdd ) * 100 else
    percent_change <- log_delta_required * 100

  return(percent_change)
}



