# find_n_after_2samp

Find the minimum number of after measurements per site needed to achieve
target power for a two-sample (unpaired) t-test that pools all before
and after measurements, ignoring within-site correlation.

## Usage

``` r
find_n_after_2samp(
  S,
  nB,
  delta,
  sd_w,
  target_power = 0.8,
  alpha = 0.05,
  n_grid = 1:50
)
```

## Arguments

- S:

  Number of sites

- nB:

  Number of before measurements per site

- delta:

  Hypothesized mean change

- sd_w:

  Within-site standard deviation

- target_power:

  Desired power level (default 0.8)

- alpha:

  Significance level (default 0.05)

- n_grid:

  Grid of after measurements to evaluate (default 1:50)

## Value

A list with `n_star` (minimum after measurements for target power) and
`curve` (data frame of n_after and power)
