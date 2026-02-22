# find_min_detectable_percent

Calculates the smallest percentage change detectable for a given sample
size.

## Usage

``` r
find_min_detectable_percent(
  S,
  nB,
  nA,
  sd_within = NA,
  sd_delta = 0,
  logTransform = FALSE,
  logAdd = 0,
  baseline_mean = NULL,
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

- sd_within:

  Within-site standard deviation (log calculated if logTransform=TRUE)

- sd_delta:

  Between-site SD of true changes (default 0)

- logTransform:

  Logical indicating whether to calculate the detectable change on the
  log scale (default FALSE)

- logAdd:

  Value to add to response variable before log-transforming to avoid
  issues

- baseline_mean:

  Mean of the original response variable before the change. Required
  when logTransform=TRUE to convert percent change back to original
  scale.

- target_power:

  Desired power (default 0.8)

- alpha:

  Significance level (default 0.05)

## Value

The minimum detectable percentage change (e.g., 30 for a 30% change)
