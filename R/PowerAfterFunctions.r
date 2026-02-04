
#' power_for_n_after
#'
#' Estimate power for a paired t-test on the average of a given number of
#' measurements at a fixed number of sites before and after a changepoint.
#'
#' @param S Number of sites
#' @param mB  Number of before measurements per site
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
power_for_n_after <- function(S, mB, nA,
                              delta, sd_w, sd_d = 0,
                              alpha = 0.05, nsim = 2000, seed = 1) {
  set.seed(seed)
  pvals <- replicate(nsim, {
    true_change <- rnorm(S, mean = delta, sd = sd_d)
    yB <- matrix(rnorm(S * mB, mean = 0, sd = sd_w), nrow = S)
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
#' @param mB Number of before measurements per site
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
find_n_after <- function(S, mB,
                         delta, sd_w, sd_d = 0,
                         target_power = 0.8,
                         alpha = 0.05,
                         n_grid = 1:50,
                         nsim = 2000, seed = 1) {
  pow <- sapply(n_grid, function(nA) {
    power_for_n_after(S = S, mB = mB, nA = nA,
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
#' @param mB Number of before measurements per site
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
power_for_sites <- function(S, mB, nA,
                            delta, sd_w, sd_d = 0,
                            alpha = 0.05, nsim = 2000, seed = 1) {
  power_for_n_after(S = S, mB = mB, nA = nA,
                    delta = delta, sd_w = sd_w, sd_d = sd_d,
                    alpha = alpha, nsim = nsim, seed = seed)
}

#' find_min_sites
#'
#' Function to find the minimum number of sites needed to achieve target power from
#' a paired t test on the means across observations at sites before and after a changepoint.
#'
#' @param mB Number of before measurements per site
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
find_min_sites <- function(mB, nA,
                           delta, sd_w, sd_d = 0,
                           target_power = 0.8,
                           alpha = 0.05,
                           S_grid = 2:50,
                           nsim = 2000, seed = 1) {
  S_grid <- S_grid[S_grid >= 2]
  pow <- sapply(S_grid, function(S) {
    power_for_sites(S = S, mB = mB, nA = nA,
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

