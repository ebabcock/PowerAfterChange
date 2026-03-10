# find_min_detectable_percent_2samp

Calculates the smallest percentage change detectable for a given sample
size using a two-sample (unpaired) t-test that ignores within-site
correlation.

## Usage

``` r
find_min_detectable_percent_2samp(
  S,
  nB,
  nA,
  sd_pooled,
  baseline_mean,
  target_power = 0.8,
  alpha = 0.05
)
```

## Arguments

- S:

  Number of sites

- nB:

  Number of before measurements per site

- nA:

  Number of after measurements per site

- sd_pooled:

  Standard deviation among all data points

- baseline_mean:

  Mean of the response variable before the change (used to convert the
  detectable absolute delta to a percent change)

- target_power:

  Desired power (default 0.8)

- alpha:

  Significance level (default 0.05)

## Value

The minimum detectable percentage change (e.g., 30 for a 30% change)
