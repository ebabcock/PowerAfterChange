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
#' @param distribution Distribution for simulated data. One of
#'   "normal", "nbinom", or "binomial".
#' @param useTest Which test to use. One of "paired-t", "wilcoxon", or "prop.test".
#' @param nbinom_mu Mean parameter for negative binomial (mu)
#' @param nbinom_disp Dispersion (size) parameter for negative binomial
#' @param binomial_size Size parameter (trials) for binomial
#' @param binomial_prob Probability parameter for binomial
#'
#' @returns Estimated power for the given number of after measurements
#' @export
#'
power_for_n_after <- function(S, nB, nA,
                              delta, sd_w, sd_d = 0,
                              alpha = 0.05, nsim = 2000, seed = 1,
                              distribution = c("normal", "nbinom", "binomial"),
                              useTest = c("paired-t", "wilcoxon", "prop.test"),
                              nbinom_mu = NULL, nbinom_disp = NULL,
                              binomial_size = NULL, binomial_prob = NULL) {
  set.seed(seed)
  distribution <- match.arg(distribution)
  useTest <- match.arg(useTest)

  if (useTest == "prop.test" && distribution != "binomial") {
    stop("useTest = 'prop.test' is only supported for distribution = 'binomial'.")
  }

  if (distribution == "nbinom" && (is.null(nbinom_mu) || is.null(nbinom_disp))) {
    stop("For distribution = 'nbinom', provide nbinom_mu and nbinom_disp.")
  }

  if (distribution == "binomial" && (is.null(binomial_size) || is.null(binomial_prob))) {
    stop("For distribution = 'binomial', provide binomial_size and binomial_prob.")
  }

  pvals <- replicate(nsim, {
    true_change <- rnorm(S, mean = delta, sd = sd_d)

    if (distribution == "normal") {
      yB <- matrix(rnorm(S * nB, mean = 0, sd = sd_w), nrow = S)
      yA <- matrix(rnorm(S * nA, mean = true_change, sd = sd_w), nrow = S)
    }

    if (distribution == "nbinom") {
      etaB <- matrix(rnorm(S * nB, mean = log(nbinom_mu), sd = sd_w), nrow = S)
      etaA <- matrix(rnorm(S * nA, mean = log(nbinom_mu) + true_change, sd = sd_w), nrow = S)
      muB <- exp(etaB)
      muA <- exp(etaA)
      yB <- matrix(rnbinom(S * nB, size = nbinom_disp, mu = muB), nrow = S)
      yA <- matrix(rnbinom(S * nA, size = nbinom_disp, mu = muA), nrow = S)
    }

    if (distribution == "binomial") {
      etaB <- matrix(rnorm(S * nB, mean = qlogis(binomial_prob), sd = sd_w), nrow = S)
      etaA <- matrix(rnorm(S * nA, mean = qlogis(binomial_prob) + true_change, sd = sd_w), nrow = S)
      probB <- plogis(etaB)
      probA <- plogis(etaA)
      yB <- matrix(rbinom(S * nB, size = binomial_size, prob = probB), nrow = S)
      yA <- matrix(rbinom(S * nA, size = binomial_size, prob = probA), nrow = S)
    }

    if (useTest %in% c("paired-t", "wilcoxon")) {
      if (distribution == "binomial") {
        D <- rowMeans(yA / binomial_size) - rowMeans(yB / binomial_size)
      } else {
        D <- rowMeans(yA) - rowMeans(yB)
      }

      if (useTest == "paired-t") {
        t.test(D, mu = 0)$p.value
      } else {
        wilcox.test(D, mu = 0)$p.value
      }
    } else {
      total_success <- c(sum(yB), sum(yA))
      total_trials <- c(S * nB * binomial_size, S * nA * binomial_size)
      prop.test(total_success, total_trials, correct = FALSE)$p.value
    }
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
#' @param distribution Distribution for simulated data. One of
#'   "normal", "nbinom", or "binomial".
#' @param useTest Which test to use. One of "paired-t", "wilcoxon", or "prop.test".
#' @param nbinom_mu Mean parameter for negative binomial (mu)
#' @param nbinom_disp Dispersion (size) parameter for negative binomial
#' @param binomial_size Size parameter (trials) for binomial
#' @param binomial_prob Probability parameter for binomial
#'
#' @returns A list containing the minimum number of after measurements needed to achieve the target power and a data frame of the power curve
#' @export
#'
find_n_after <- function(S, nB,
                         delta, sd_w, sd_d = 0,
                         target_power = 0.8,
                         alpha = 0.05,
                         n_grid = 1:50,
                         nsim = 2000, seed = 1,
                         distribution = c("normal", "nbinom", "binomial"),
                         useTest = c("paired-t", "wilcoxon", "prop.test"),
                         nbinom_mu = NULL, nbinom_disp = NULL,
                         binomial_size = NULL, binomial_prob = NULL) {
  pow <- sapply(n_grid, function(nA) {
    power_for_n_after(S = S, nB = nB, nA = nA,
                      delta = delta, sd_w = sd_w, sd_d = sd_d,
                      alpha = alpha, nsim = nsim, seed = seed,
                      distribution = distribution,
                      useTest = useTest,
                      nbinom_mu = nbinom_mu,
                      nbinom_disp = nbinom_disp,
                      binomial_size = binomial_size,
                      binomial_prob = binomial_prob)
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
#' @param distribution Distribution for simulated data. One of
#'   "normal", "nbinom", or "binomial".
#' @param useTest Which test to use. One of "paired-t", "wilcoxon", or "prop.test".
#' @param nbinom_mu Mean parameter for negative binomial (mu)
#' @param nbinom_disp Dispersion (size) parameter for negative binomial
#' @param binomial_size Size parameter (trials) for binomial
#' @param binomial_prob Probability parameter for binomial
#'
#' @returns Estimated power for the given number of sites
#' @export
#'
power_for_sites <- function(S, nB, nA,
                            delta, sd_w, sd_d = 0,
                            alpha = 0.05, nsim = 2000, seed = 1,
                            distribution = c("normal", "nbinom", "binomial"),
                            useTest = c("paired-t", "wilcoxon", "prop.test"),
                            nbinom_mu = NULL, nbinom_disp = NULL,
                            binomial_size = NULL, binomial_prob = NULL) {
  power_for_n_after(S = S, nB = nB, nA = nA,
                    delta = delta, sd_w = sd_w, sd_d = sd_d,
                    alpha = alpha, nsim = nsim, seed = seed,
                    distribution = distribution,
                    useTest = useTest,
                    nbinom_mu = nbinom_mu,
                    nbinom_disp = nbinom_disp,
                    binomial_size = binomial_size,
                    binomial_prob = binomial_prob)
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
#' @param distribution Distribution for simulated data. One of
#'   "normal", "nbinom", or "binomial".
#' @param useTest Which test to use. One of "paired-t", "wilcoxon", or "prop.test".
#' @param nbinom_mu Mean parameter for negative binomial (mu)
#' @param nbinom_disp Dispersion (size) parameter for negative binomial
#' @param binomial_size Size parameter (trials) for binomial
#' @param binomial_prob Probability parameter for binomial
#'
#' @returns A list containing the minimum number of sites needed to achieve the target power and a data frame of the power curve
#' @export
#'
find_min_sites <- function(nB, nA,
                           delta, sd_w, sd_d = 0,
                           target_power = 0.8,
                           alpha = 0.05,
                           S_grid = 2:50,
                           nsim = 2000, seed = 1,
                           distribution = c("normal", "nbinom", "binomial"),
                           useTest = c("paired-t", "wilcoxon", "prop.test"),
                           nbinom_mu = NULL, nbinom_disp = NULL,
                           binomial_size = NULL, binomial_prob = NULL) {
  S_grid <- S_grid[S_grid >= 2]
  pow <- sapply(S_grid, function(S) {
    power_for_sites(S = S, nB = nB, nA = nA,
                    delta = delta, sd_w = sd_w, sd_d = sd_d,
                    alpha = alpha, nsim = nsim, seed = seed,
                    distribution = distribution,
                    useTest = useTest,
                    nbinom_mu = nbinom_mu,
                    nbinom_disp = nbinom_disp,
                    binomial_size = binomial_size,
                    binomial_prob = binomial_prob)
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
#' @param S_grid Number of sites grid
#' @param delta_target Target mean change
#' @param sd_diff Standard deviation of difference in means
#'
#' @returns Power for each number of sites in S_grid
#' @export
find_min_sites_analytical <- function(S_grid, delta_target, sd_diff) {
  S_grid <- S_grid[S_grid >= 2]
  sapply(S_grid, function(n) {
    power.t.test(
      n = n,
      delta = delta_target,
      sd = sd_diff,
      sig.level = 0.05,
      type = "paired",
      alternative = "two.sided"
    )$power
  })
}


#' power_for_nA
#'
#' Function to calculate power for given number of after samples for a paired t-test
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
power_for_nA_analytical <- function(nA, S, nB, delta, sd_within, sd_delta, alpha) {
  sd_diff <- getSD_difference(sd_within, nA, nB, sd_delta)
  power_result <- power.t.test(
    n = S,
    delta = delta,
    sd = sd_diff,
    sig.level = alpha,
    power = NULL,
    type = "paired",
    alternative = "two.sided"
  )
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
#' @param groupVar Optional character vector of grouping variables to summarize by
#' @param typeTansform Character indicating transformation the response variable, "none","log", or "sqrt" (default "none")"
#' @param addValue Value to add to response variable before transforming (default 0)
#'
#' @returns A data frame containing the estimated mean and standard deviation of the
#' response variable across all sites and data points (grand_mean,grand_sd), and the
#' within and across site standard deviations (sd_within,
#' sd_between) in the original scale, and optionally the log or sqrt transformed
#' scale if typeTransform is "log" or "sqrt", with an added constant if
#' provide d(addValue). The number of sites in the baseline is provided (nB), or, if
#' the number of samples per site is not consistent, the number of complete
#' and the number of sites with fewer samples (nb_complete and nB_incomplete).
#' The proportion of positive values in the response variable is also calculated (prop_positive) to help assess whether a log transformation is appropriate.
#' If groupVar is provided, these statistics are calculated separately for each group defined by groupVar.
#'
#' @export
#'
summarize_baseline <- function(baseline,
                               siteVar = "site",
                               responseVar = "y",
                               groupVar = NULL,
                               typeTransform = "none",
                               addValue = 0) {
  group_syms <- if (is.null(groupVar)) character(0) else groupVar
  if(typeTransform == "log"){
    logTransform <- TRUE
    sqrtTransform <- FALSE
  } else if (typeTransform == "sqrt"){
    logTransform <- FALSE
    sqrtTransform <- TRUE
  } else {
    logTransform <- FALSE
    sqrtTransform <- FALSE
  }
  returnVal <- baseline %>%
    group_by(across(all_of(c(group_syms)))) %>%
    mutate(grand_mean=mean(!!sym(responseVar)),
           grand_sd=sd(!!sym(responseVar)),
           grand_logmean=if_else(logTransform, mean(log(!!sym(responseVar) + addValue)), NA_real_),
           grand_logsd=if_else(logTransform, sd(log(!!sym(responseVar) + addValue)), NA_real_),
           grand_sqrtmean=if_else(sqrtTransform, mean(sqrt(!!sym(responseVar) + addValue)), NA_real_),
           grand_sqrtsd=if_else(sqrtTransform, sd(sqrt(!!sym(responseVar) + addValue)), NA_real_)
           ) %>%
    group_by(across(all_of(c(group_syms, siteVar)))) %>%
    summarize(
      grand_mean = unique(grand_mean),
      grand_sd = unique(grand_sd),
      grand_logmean = unique(grand_logmean),
      grand_logsd = unique(grand_logsd),
      grand_sqrtmean = unique(grand_sqrtmean),
      grand_sqrtsd = unique(grand_sqrtsd),
      site_mean = mean(!!sym(responseVar)),
      site_sd = sd(!!sym(responseVar)),
      site_logmean = if_else(logTransform, mean(log(!!sym(responseVar) + addValue)), NA_real_),
      site_logsd = if_else(logTransform, sd(log(!!sym(responseVar) + addValue)), NA_real_),
      site_sqrtmean = if_else(sqrtTransform, mean(sqrt(!!sym(responseVar) + addValue)), NA_real_),
      site_sqrtsd = if_else(sqrtTransform, sd(sqrt(!!sym(responseVar) + addValue)), NA_real_),
      prop_positive_site = mean(!!sym(responseVar) > 0),
      n = n(),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(group_syms))) %>%
    summarize(
      grand_mean = unique(grand_mean),
      grand_sd = unique(grand_sd),
      grand_logmean = unique(grand_logmean),
      grand_logsd = unique(grand_logsd),
      grand_sqrtmean = unique(grand_sqrtmean),
      grand_sqrtsd = unique(grand_sqrtsd),
      sd_within = mean(site_sd),
      sd_between = sd(site_mean),
      logsd_within = if_else(logTransform, mean(site_logsd), NA_real_),
      logsd_between = if_else(logTransform, sd(site_logmean), NA_real_),
      sqrtsd_within = if_else(sqrtTransform, mean(site_sqrtsd), NA_real_),
      sqrtsd_between = if_else(sqrtTransform, sd(site_sqrtmean), NA_real_),
      prop_positive = mean(prop_positive_site),
      n_sites = n(),
      nB = if_else(length(unique(n)) == 1, unique(n)[1], NA_real_),
      nB_complete= if_else(length(unique(n)) >1 ,min(n), NA_real_),
      nB_incomplete = if_else(length(unique(n))>1 ,max(n) - min(n), NA_real_),
      .groups = "drop"
    )%>%
    select(where(~ any(!is.na(.))))

  if (any(is.na(returnVal$nB))) {
    warning("Number of before samples per site is not consistent across sites. nB is set to NA.")
  }
  return(returnVal)
}


#' find_desired_change
#'
#' Function to convert a desired change in the original response variable, specified
#' as either an absolute change or a percent change, to the corresponding change
#' on the scale of the response variable (e.g. log or sqrt) if a transformation is
#' specified, using the baseline mean to convert percent changes to absolute
#' changes on the original scale before applying the transformation.
#'
#' @param change_type Character indicating whether the desired change is specified
#' as an "absolute" change or a "percent" change in the original response variable.
#' @param change_value Numeric value of the desired change in the original response.
#' For example, if change_type is "percent", a value of 30 means a desired change
#' of 30% (i.e. multiply original mean by 1.3).
#' @param baseline_mean Mean of the original response variable before the change
#' @param typeTransform Character indicating the transformation to apply to the response
#'  variable, one of "none", "log", or "sqrt" (default "none").
#' @param addValue Value to add to response variable before transforming (default 0).
#'
#' @export
#'
find_desired_change <- function(change_type,
                                change_value,
                                baseline_mean,
                                typeTransform = c("none", "log", "sqrt"),
                                addValue = 0) {
  typeTransform <- match.arg(typeTransform)
  # convert percent change to absolute change on original scale
  if (change_type == "percent") {
    abs_change <- baseline_mean * (change_value / 100)
  } else if (change_type == "absolute") {
    abs_change <- change_value
  } else {
    stop("change_type must be either 'absolute' or 'percent'.")
  }
  # baseline and changed value on original scale
  baseline <- baseline_mean
  changed <- baseline_mean + abs_change
  # apply transformation (with addValue shift)
  transform_fun <- switch(
    typeTransform,
    none = function(x) x,
    log  = function(x) log(x),
    sqrt = function(x) sqrt(x)
  )
  baseline_t <- transform_fun(baseline + addValue)
  changed_t  <- transform_fun(changed + addValue)
  # change on transformed scale
  changed_t - baseline_t
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
#' @param typeTransform Character indicating the transformation to apply to the response
#'  variable, one of "none", "log", or "sqrt" (default "none").
#' @param addValue Value to add to response variable before transforming to avoid issues
#' @param baseline_mean Mean of the original response variable before the change.
#'   Required when typeTransform is "log" or "sqrt" to convert percent change back to original scale.
#' @param target_power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#'
#' @returns The minimum detectable percentage change (e.g., 30 for a 30% change)
#' @export
find_min_detectable_percent <- function(S,
                                        nB,
                                        nA,
                                        sd_within = NA,
                                        sd_delta = 0,
                                        typeTransform = c("none", "log", "sqrt"),
                                        addValue = 0,
                                        baseline_mean = NULL,
                                        target_power = 0.8,
                                        alpha = 0.05) {

  typeTransform <- match.arg(typeTransform)

  if (typeTransform %in% c("log", "sqrt") && is.null(baseline_mean)) {
    stop("baseline_mean must be provided when typeTransform is 'log' or 'sqrt'.")
  }
  if (typeTransform == "log" && baseline_mean + addValue <= 0) {
    stop("baseline_mean + addValue must be > 0 when typeTransform is 'log'.")
  }
  if (typeTransform == "sqrt" && baseline_mean + addValue < 0) {
    stop("baseline_mean + addValue must be >= 0 when typeTransform is 'sqrt'.")
  }

  power_root_func <- function(delta) {
    sd_diff <- getSD_difference(sd_within, nA, nB, sd_delta)

    current_power <- power.t.test(
      n = S,
      delta = delta,
      sd = sd_diff,
      sig.level = alpha,
      type = "paired",
      alternative = "two.sided"
    )$power

    current_power - target_power
  }

  fit <- uniroot(power_root_func, interval = c(1e-4, 5))
  delta_required <- fit$root

  if (typeTransform == "none") {
    percent_change <- delta_required * 100
  } else {
    transform_fun <- switch(
      typeTransform,
      log  = function(x) log(x),
      sqrt = function(x) sqrt(x)
    )
    inv_transform_fun <- switch(
      typeTransform,
      log  = function(x) exp(x),
      sqrt = function(x) x^2
    )

    baseline_t <- transform_fun(baseline_mean + addValue)
    changed_t  <- baseline_t + delta_required

    before <- baseline_mean
    after  <- inv_transform_fun(changed_t) - addValue

    percent_change <- ((after - before) / before) * 100
  }

  return(percent_change)
}

#' power_for_percent_change
#'
#' Function to calculate power for a given percent change in the mean
#' before vs. after change using power.t.test (paired t-test on site-level differences).
#'
#' @param percent_change Target percent change in the original response variable
#' @param baseline_mean Mean of the original response variable before the change
#' @param nA Number of after samples per site
#' @param S Number of sites
#' @param nB Number of before samples per site
#' @param sd_within Standard deviation within sites (in transformed scale)
#' @param sd_delta Standard deviation of true changes among sites
#' @param alpha Significance level
#' @param typeTransform Character indicating the transformation to apply to the response
#'  variable, one of "none", "log", or "sqrt" (default "none").
#' @param addValue Value to add to response variable before transforming to avoid issues
#'
#' @returns Power for the given percent change
#' @export
#'
power_for_percent_change <- function(percent_change,
                                     baseline_mean = NULL,
                                     nA,
                                     S,
                                     nB,
                                     sd_within,
                                     sd_delta,
                                     alpha,
                                     typeTransform = c("none", "log", "sqrt"),
                                     addValue = 0) {
  typeTransform <- match.arg(typeTransform)

  if (is.null(baseline_mean)) {
    stop("baseline_mean must be provided.")
  }
  if (typeTransform == "log" && baseline_mean + addValue <= 0) {
    stop("baseline_mean + addValue must be > 0 when typeTransform is 'log'.")
  }
  if (typeTransform == "sqrt" && baseline_mean + addValue < 0) {
    stop("baseline_mean + addValue must be >= 0 when typeTransform is 'sqrt'.")
  }

  before <- baseline_mean
  after <- baseline_mean * (1 + percent_change / 100)

  if (typeTransform == "none") {
    delta <- after - before
  } else {
    transform_fun <- switch(
      typeTransform,
      log  = function(x) log(x),
      sqrt = function(x) sqrt(x)
    )
    delta <- transform_fun(after + addValue) - transform_fun(before + addValue)
  }

  sd_diff <- getSD_difference(sd_within, nA, nB, sd_delta)

  power_result <- power.t.test(
    n = S,
    delta = delta,
    sd = sd_diff,
    sig.level = alpha,
    power = NULL,
    type = "paired",
    alternative = "two.sided"
  )

  return(power_result$power)
}

#' power_2samp_analytical
#'
#' Analytical power for a two-sample (unpaired) t-test
#' using the non-central t distribution. Handles unbalanced group sizes.
#'
#' @param n_before Number of before measurements (total across all sites)
#' @param n_after Number of after measurements (total across all sites)
#' @param delta Hypothesized mean change (absolute, on the scale of the response)
#' @param sd_w Within-site standard deviation (or pooled standard deviation across all samples)
#' @param alpha Significance level (default 0.05)
#'
#' @export
#'
power_2samp_analytical <- function(n_before, n_after, delta, sd_pooled, alpha) {
  df    <- n_before + n_after - 2
  ncp   <- delta / (sd_pooled * sqrt(1 / n_before + 1 / n_after))
  t_crit <- qt(1 - alpha / 2, df = df)
  pt(-t_crit, df = df, ncp = ncp) +
    pt(t_crit, df = df, ncp = ncp, lower.tail = FALSE)
}

#' power_for_n_after_2samp
#'
#' Estimate power for a two-sample (unpaired) t-test comparing all before vs. all
#' after measurements, ignoring possible within-site correlation.
#'
#' @param S Number of sites
#' @param nB Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param delta Hypothesized mean change
#' @param sd_pooled Standard deviation between samples
#' @param sd_d Between-site standard deviation of true changes (default 0)
#' @param alpha Significance level, defaults to 0.05
#' @param nsim Number of simulations to run
#' @param seed Random seed for reproducibility
#'
#' @returns Estimated power (numeric scalar)
#' @export
power_for_n_after_2samp <- function(S, nB, nA,
                                    delta, sd_pooled, sd_d = 0,
                                    alpha = 0.05, nsim = 2000, seed = 1) {
  set.seed(seed)
  pvals <- replicate(nsim, {
    true_change <- rnorm(S, mean = delta, sd = sd_d)
    yB <- matrix(rnorm(S * nB, mean = 0, sd = sd_pooled), nrow = S)
    yA <- matrix(rnorm(S * nA, mean = true_change, sd = sd_pooled), nrow = S)
    all_before <- as.vector(yB)
    all_after  <- as.vector(yA)
    t.test(all_after, all_before, var.equal = TRUE)$p.value
  })
  mean(pvals < alpha)
}

#' find_min_detectable_percent_2samp
#'
#' Calculates the smallest percentage change detectable for a given sample size
#' using a two-sample (unpaired) t-test that ignores within-site correlation.
#'
#' @param S Number of sites
#' @param nB Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param sd_pooled Standard deviation among all data points
#' @param baseline_mean Mean of the response variable before the change (used to
#'   convert the detectable absolute delta to a percent change)
#' @param target_power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#'
#' @returns The minimum detectable percentage change (e.g., 30 for a 30% change)
#' @export
find_min_detectable_percent_2samp <- function(S, nB, nA,
                                              sd_pooled,
                                              baseline_mean,
                                              target_power = 0.8,
                                              alpha = 0.05) {
  n_before <- S * nB
  n_after  <- S * nA
  power_root_func <- function(delta) {
    power_2samp_analytical(n_before, n_after, delta, sd_pooled, alpha) - target_power
  }
  fit <- uniroot(power_root_func, interval = c(1e-6, 1e6))  # wide enough for any realistic delta
  delta_required <- fit$root
  (delta_required / baseline_mean) * 100
}

#' find_min_sites_2samp
#'
#' Find the minimum number of sites needed to achieve target power for a
#' two-sample (unpaired) t-test that pools all before and after measurements,
#' ignoring within-site correlation.
#'
#' @param nB Number of before measurements per site
#' @param nA Number of after measurements per site
#' @param delta Hypothesized mean change (absolute, on the scale of the response)
#' @param sd_pooled Standard deviation among all samples
#' @param target_power Desired power level (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param S_grid Grid of site numbers to evaluate (default 2:50)
#'
#' @returns A list with `S_star` (minimum sites for target power) and `curve`
#'   (data frame of S and power)
#' @export
#'
find_min_sites_2samp <- function(nB, nA,
                                 delta, sd_pooled,
                                 target_power = 0.8,
                                 alpha = 0.05,
                                 S_grid = 2:50) {
  pow <- sapply(S_grid, function(S) {
    power_2samp_analytical(S * nB, S * nA, delta, sd_pooled, alpha)
  })
  out    <- data.frame(S = S_grid, power = pow)
  S_star <- out$S[which(out$power >= target_power)[1]]
  list(S_star = S_star, curve = out)
}

#' find_n_after_2samp
#'
#' Find the minimum number of after measurements per site needed to achieve
#' target power for a two-sample (unpaired) t-test that pools all before and
#' after measurements, ignoring within-site correlation.
#'
#' @param S Number of sites
#' @param nB Number of before measurements per site
#' @param delta Hypothesized mean change
#' @param sd_pooled Standard deviation among all samples
#' @param target_power Desired power level (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param n_grid Grid of after measurements to evaluate (default 1:50)
#'
#' @returns A list with `n_star` (minimum after measurements for target power)
#'   and `curve` (data frame of n_after and power)
#' @export
find_n_after_2samp <- function(S, nB,
                               delta,
                               sd_pooled,
                               target_power = 0.8,
                               alpha = 0.05,
                               n_grid = 1:50) {
  pow <- sapply(n_grid, function(nA) {
    power_2samp_analytical(S * nB, S * nA, delta, sd_pooled, alpha)
  })
  out    <- data.frame(n_after = n_grid, power = pow)
  n_star <- out$n_after[which(out$power >= target_power)[1]]
  list(n_star = n_star, curve = out)
}
