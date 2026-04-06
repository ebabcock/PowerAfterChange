# find_min_sites_2samp

Find the minimum number of sites needed to achieve target power for a
two-sample (unpaired) t-test that pools all before and after
measurements, ignoring within-site correlation.

## Usage

``` r
find_min_sites_2samp(
  nB,
  nA,
  n_sites_before = NULL,
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

- S_before:

  Number of sites before, if you wish to keep this number constant in
  the analysis. Defaults to NULL, which allows the number of sites
  before and after to vary in the calculation. If n_sites_before is
  provided, the function will calculate power for a fixed number of
  sites before and varying number of sites after.

## Value

A list with `S_star` (minimum sites for target power) and `curve` (data
frame of S and power). If nB and nA are the number of years before
after, then S_star is the minimum number of sites needed to achieve
target power with a sample size of `n_before=nB*S` and `n_before=nA*S`
in a 2-sample t-test. If n_sites_before is provided, then S_star is the
minimum number of sites needed to achieve target power with S_before
sites before and S_star sites after.
