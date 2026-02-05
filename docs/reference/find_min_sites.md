# find_min_sites

Function to find the minimum number of sites needed to achieve target
power from a paired t test on the means across observations at sites
before and after a changepoint.

## Usage

``` r
find_min_sites(
  mB,
  nA,
  delta,
  sd_w,
  sd_d = 0,
  target_power = 0.8,
  alpha = 0.05,
  S_grid = 2:50,
  nsim = 2000,
  seed = 1
)
```

## Arguments

- nA:

  Number of after measurements per site

- delta:

  Hypothesized mean change

- sd_w:

  Within-site standard deviation

- sd_d:

  Between-site standard deviation of true changes

- target_power:

  Desired power level to achieve

- alpha:

  Significance level, defaults to 0.05

- S_grid:

  Grid of site numbers to evaluate

- nsim:

  Number of simulations to run

- seed:

  Random seed for reproducibility

- nB:

  Number of before measurements per site

## Value

A list containing the minimum number of sites needed to achieve the
target power and a data frame of the power curve
