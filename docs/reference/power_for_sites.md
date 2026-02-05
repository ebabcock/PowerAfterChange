# power_for_sites

Function to estimate power for a paired t-test on the average of a given
number of sites before and after a changepoint.

## Usage

``` r
power_for_sites(
  S,
  mB,
  nA,
  delta,
  sd_w,
  sd_d = 0,
  alpha = 0.05,
  nsim = 2000,
  seed = 1
)
```

## Arguments

- S:

  Number of sites

- mB:

  Number of before measurements per site

- nA:

  Number of after measurements per site

- delta:

  Hypothesized mean change

- sd_w:

  Within-site standard deviation

- sd_d:

  Between-site standard deviation of true changes

- alpha:

  Significance level, defaults to 0.05

- nsim:

  Number of simulations to run

- seed:

  Random seed for reproducibility

## Value

Estimated power for the given number of sites
