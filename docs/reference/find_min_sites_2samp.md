# find_min_sites_2samp

Find the minimum number of sites needed to achieve target power for a
two-sample (unpaired) t-test that pools all before and after
measurements, ignoring within-site correlation.

## Usage

``` r
find_min_sites_2samp(
  nB,
  nA,
  delta,
  sd_pooled,
  target_power = 0.8,
  alpha = 0.05,
  S_grid = 2:50
)
```

## Arguments

- nB:

  Number of before measurements per site

- nA:

  Number of after measurements per site

- delta:

  Hypothesized mean change (absolute, on the scale of the response)

- sd_pooled:

  Standard deviation among all samples

- target_power:

  Desired power level (default 0.8)

- alpha:

  Significance level (default 0.05)

- S_grid:

  Grid of site numbers to evaluate (default 2:50)

## Value

A list with `S_star` (minimum sites for target power) and `curve` (data
frame of S and power)
