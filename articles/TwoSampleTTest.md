# Two-sample t-test power analysis

## Overview

This is a work in progress, the functions are not final.

The functions in this vignette calculate power and sample size for a
**two-sample (unpaired) t-test** that compares all before observations
pooled together against all after observations pooled together, without
accounting for any repeated-measures structure of the data (i.e.,
multiple samples taken at the same site before and after the change).
This would be correct if the data points were collected at random sites.
It would be approximately valid if there was not much repeatability in
the measurments at each site in a repeated measures design.

The functions are:

- [`power_for_n_after_2samp()`](https://ebabcock.github.io/PowerAfterChange/reference/power_for_n_after_2samp.md)
  — simulation-based power for a two-sample t-test
- [`find_min_sites_2samp()`](https://ebabcock.github.io/PowerAfterChange/reference/find_min_sites_2samp.md)
  — analytical minimum sites needed for target power
- [`find_n_after_2samp()`](https://ebabcock.github.io/PowerAfterChange/reference/find_n_after_2samp.md)
  — analytical minimum after samples per site for target power
- [`find_min_detectable_percent_2samp()`](https://ebabcock.github.io/PowerAfterChange/reference/find_min_detectable_percent_2samp.md)
  — analytical minimum detectable percent change

``` r
library(tidyverse)
theme_set(theme_minimal())
library(PowerAfterChange)
```

## Simulated data

We use the **same simulated dataset** as `UserGuide.Rmd` so that results
can be directly compared.

``` r
set.seed(123)
S_demo    <- 12
nB_demo   <- 5
sd_within <- 2
sd_between <- 1
before_mean <- 10
siteMean <- rnorm(S_demo, mean = before_mean, sd = sd_between)
baseline_demo <- data.frame(
  site = rep(1:S_demo, each = nB_demo)) %>%
  mutate(y = rnorm(S_demo * nB_demo, mean = siteMean[site], sd = sd_within))
```

### Estimate total standard deviation among samples

``` r
sd_total <- sd(baseline_demo$y)
```

### Planning assumptions

``` r
delta_target <- 1    # absolute change to detect
sd_delta     <- 0.5  # site-to-site variability in true changes
```

------------------------------------------------------------------------

## Question 1: Minimum number of sites (nB = 5, nA = 5)

``` r
sites_2samp <- find_min_sites_2samp(
  nB = nB_demo, nA = 5,
  delta = delta_target,
  sd_pooled  = sd_total,
  target_power = 0.8, alpha = 0.05,
  S_grid = 2:50
)
sites_2samp$S_star
```

    ## [1] 12

``` r
ggplot(sites_2samp$curve, aes(x = S, y = power)) +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_vline(xintercept = sites_2samp$S_star, linetype = "dashed", color = "blue") +
  labs(title = "Two-sample t-test: Power Curve for Number of Sites",
       x = "Number of Sites",
       y = "Power")
```

![](TwoSampleTTest_files/figure-html/unnamed-chunk-6-1.png)

## Question 2: Minimum after samples per site (S = 12, nB = 5)

``` r
n_after_2samp <- find_n_after_2samp(
  S  = S_demo, nB = nB_demo,
  delta = delta_target,
  sd_pooled  = sd_total,
  target_power = 0.8, alpha = 0.05,
  n_grid = 1:40
)
n_after_2samp$n_star
```

    ## [1] 5

``` r
ggplot(n_after_2samp$curve, aes(x = n_after, y = power)) +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_vline(xintercept = n_after_2samp$n_star, linetype = "dashed", color = "blue") +
  labs(title = "Two-sample t-test: Power Curve for After Samples per Site",
       x = "Number of After Samples per Site",
       y = "Power")
```

![](TwoSampleTTest_files/figure-html/unnamed-chunk-8-1.png)

------------------------------------------------------------------------

## Question 3: Minimum detectable percent change (S = 12, nB = 5, nA = 5)

``` r
min_pct_2samp <- find_min_detectable_percent_2samp(
  S = S_demo, nB = nB_demo, nA = 5,
  sd_pooled = sd_total,
  baseline_mean = before_mean,
  target_power = 0.8, alpha = 0.05
)
min_pct_2samp
```

    ## [1] 9.813063

How does the minimum detectable percent change vary with the number of
after samples?

``` r
nA_grid <- 1:40
detectable_2samp_df <- data.frame(
  n_after = nA_grid,
  min_detectable_percent = sapply(nA_grid, function(nA) {
    find_min_detectable_percent_2samp(
      S = S_demo, nB = nB_demo, nA = nA,
      sd_pooled = sd_total,
      baseline_mean = before_mean,
      target_power = 0.8, alpha = 0.05
    )
  })
)

ggplot(detectable_2samp_df, aes(x = n_after, y = min_detectable_percent)) +
  geom_line() +
  geom_hline(yintercept = min_pct_2samp, linetype = "dashed", color = "red") +
  labs(
    title = "Two-sample t-test: Minimum Detectable Percent Change vs. After Samples",
    x = "Number of After Samples per Site",
    y = "Minimum Detectable Percent Change"
  )
```

![](TwoSampleTTest_files/figure-html/unnamed-chunk-10-1.png)

------------------------------------------------------------------------

## Comparing paired t-test vs. two-sample t-test power curves

The plot below overlays the power curve from the **paired t-test**
(`power_for_nA_analytical`) against the **two-sample t-test**
(`find_n_after_2samp`) for the same design parameters. The paired design
accounts for between-site variability in changes (`sd_d`), while the
two-sample design ignores site structure entirely.

``` r
nA_grid <- 1:40

# Paired analytical power (includes sd_delta for between-site change variability)
sd_within_hat <-getSD_within(baseline = baseline_demo,
                             siteVar = "site",
                             responseVar = "y")
power_paired <- sapply(nA_grid, function(nA) {
  power_for_nA_analytical(nA, S = S_demo, nB = nB_demo,
                          delta = delta_target,
                          sd_within = sd_within_hat,
                          sd_delta = sd_delta,
                          alpha = 0.05)
})

# Two-sample analytical power (ignores site structure)
power_2samp <- n_after_2samp$curve$power

comparison_df <- data.frame(
  n_after      = nA_grid,
  Paired       = power_paired,
  Two_sample   = power_2samp
) %>%
  pivot_longer(-n_after, names_to = "Method", values_to = "power") %>%
  mutate(Method = recode(Method, Two_sample = "Two-sample"))

ggplot(comparison_df, aes(x = n_after, y = power, color = Method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
  labs(
    title = "Paired t-test vs. Two-sample t-test: Power Curves",
    subtitle = paste0("S=", S_demo, ", nB=", nB_demo,
                      ", delta=", delta_target,
                      ", sd_w=", round(sd_within_hat, 2)),
    x = "Number of After Samples per Site",
    y = "Power",
    color = "Method"
  ) +
  theme(legend.position = "bottom")
```

![](TwoSampleTTest_files/figure-html/unnamed-chunk-11-1.png)

The paired design should more efficient here: it typically reaches
target power with fewer after samples because it eliminates between-site
baseline variability as a source of error. Still working on why it
isn’t.

------------------------------------------------------------------------

## Simulation validation

We validate the analytical `find_n_after_2samp` results by comparing
them to simulation-based estimates from `power_for_n_after_2samp` for a
few values of `nA`.

``` r
nA_check <- c(2, 5, 10, 20)

sim_power <- sapply(nA_check, function(nA) {
  power_for_n_after_2samp(
    S = S_demo, nB = nB_demo, nA = nA,
    delta = delta_target, sd_pooled = sd_total, sd_d = 0,
    alpha = 0.05, nsim = 2000, seed = 42
  )
})

analytical_power <- sapply(nA_check, function(nA) {
  n_after_2samp$curve$power[n_after_2samp$curve$n_after == nA]
})

validation_df <- data.frame(
  nA              = nA_check,
  Analytical      = round(analytical_power, 3),
  Simulation      = round(sim_power, 3)
)
knitr::kable(validation_df,
             caption = "Two-sample t-test: Analytical vs. Simulation Power")
```

|  nA | Analytical | Simulation |
|----:|-----------:|-----------:|
|   2 |      0.576 |      0.561 |
|   5 |      0.815 |      0.815 |
|  10 |      0.911 |      0.912 |
|  20 |      0.952 |      0.949 |

Two-sample t-test: Analytical vs. Simulation Power

The simulation and analytical results agree closely, confirming that the
analytical formula based on the non-central t distribution is correct.

------------------------------------------------------------------------

## Conclusion

The **paired t-test** (existing functions) is preferred when sites are
genuinely paired and repeated samples are taken at each site, because it
accounts for within-site correlation and is more powerful. The
**two-sample approach** is appropriate when before and after samples
come from independent groups with no repeated site structure, or as a
conservative check or comparison. When there is positive between-site
variation in baselines, the two-sample approach will generally require
more sites or more samples to achieve the same power as the paired
design.
