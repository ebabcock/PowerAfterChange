# Using sqrt transform

Demo: Power for a 30% change using `typeTransform = "sqrt"` and
`addValue = 0.5`

``` r
library(tidyverse)
theme_set(theme_minimal())
#devtools::install_github("ebabcock/PowerAfterChange")
library(PowerAfterChange) #Library with functions for analysis
```

### Simulate some data

``` r
set.seed(123) #for reproducibility
S_demo <- 12 #Number of sites
nB_demo <- 5 #Number of before samples per site
sd_within <- 2 #within-site standard deviation
sd_between <- 1 #between-site standard deviation
before_mean <- 10 #overall mean before change
#simulate site to site variability
siteMean <- rnorm(S_demo, mean = before_mean, sd = sd_between)
#Generate baseline data (before change)
baseline_demo <- data.frame(
  site = rep(1:S_demo, each = nB_demo)) %>%
  mutate(y = rnorm(S_demo * nB_demo, mean = siteMean[site], sd = sd_within))
```

1.  Summarize baseline to get mean and within-site SD on sqrt scale

``` r
base_summary <- summarize_baseline(
  baseline=baseline_demo,
  siteVar = "site",
  responseVar = "y",
  typeTransform = "sqrt",
  addValue = 0.5
)
```

    ## Warning: There were 2 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `grand_asinmean = if_else(asinTransform, mean(asin(sqrt(y))),
    ##   NA_real_)`.
    ## Caused by warning in `asin()`:
    ## ! NaNs produced
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

    ## Warning: There were 24 warnings in `summarize()`.
    ## The first warning was:
    ## ℹ In argument: `site_asinmean = if_else(asinTransform, mean(asin(sqrt(y))),
    ##   NA_real_)`.
    ## ℹ In group 1: `site = 1`.
    ## Caused by warning in `asin()`:
    ## ! NaNs produced
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 23 remaining warnings.

Extract inputs

``` r
baseline_mean <- base_summary$grand_mean
sd_within_sqrt <- base_summary$sqrtsd_within
```

Compute power for a 30% change

``` r
power_30pct <- power_for_percent_change(
  percent_change = 30,
  baseline_mean = baseline_mean,
  nA = 5,
  S = 12,
  nB = 5,
  sd_within = sd_within_sqrt,
  sd_delta = 0,
  alpha = 0.05,
  typeTransform = "sqrt",
  addValue = 0.5
)

power_30pct
```

    ## [1] 1
