# power_for_n_after_2samp

Estimate power for a two-sample (unpaired) t-test comparing all before
vs. all after measurements, ignoring possible within-site correlation.

## Usage

``` r
power_for_n_after_2samp(
  S,
  nB,
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

- nB:

  Number of before measurements per site

- nA:

  Number of after measurements per site

- delta:

  Hypothesized mean change

- sd_w:

  Within-site standard deviation

- sd_d:

  Between-site standard deviation of true changes (default 0)

- alpha:

  Significance level, defaults to 0.05

- nsim:

  Number of simulations to run

- seed:

  Random seed for reproducibility

## Value

Estimated power (numeric scalar)
