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
  typeTransform = c("none", "log", "sqrt"),
  addValue = 0,
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

- typeTransform:

  Character indicating the transformation to apply to the response
  variable, one of "none", "log", or "sqrt" (default "none").

- addValue:

  Value to add to response variable before transforming to avoid issues

- baseline_mean:

  Mean of the original response variable before the change. Required
  when typeTransform is "log" or "sqrt" to convert percent change back
  to original scale.

- target_power:

  Desired power (default 0.8)

- alpha:

  Significance level (default 0.05)

## Value

The minimum detectable percentage change (e.g., 30 for a 30% change)
