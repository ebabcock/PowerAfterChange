# find_n_after

Function to find the minimum number of after measurements needed to
achieve target power, for a paired t-test on the before and after means.

## Usage

``` r
find_n_after(
  S,
  mB,
  delta,
  sd_w,
  sd_d = 0,
  target_power = 0.8,
  alpha = 0.05,
  n_grid = 1:50,
  nsim = 2000,
  seed = 1
)
```

## Arguments

- S:

  Number of sites

- mB:

  Number of before measurements per site

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

- n_grid:

  Grid of after measurements to evaluate

- nsim:

  Number of simulations to run

- seed:

  Random seed for reproducibility

## Value

A list containing the minimum number of after measurements needed to
achieve the target power and a data frame of the power curve
