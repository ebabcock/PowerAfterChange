# find_min_detectable_percent_2samp

Calculates the smallest percentage change detectable for a given sample
size using a two-sample (unpaired) t-test that ignores within-site
correlation.

## Usage

``` r
find_min_detectable_percent_2samp(
  S,
  S_before = NULL,
  nB,
  nA,
  sd_pooled,
  baseline_mean,
  target_power = 0.8,
  alpha = 0.05,
  typeTransform = c("none", "log", "sqrt", "asin"),
  addValue = 0
)
```

## Arguments

- S:

  Number of sites both before and after (if S_before is NULL) or number
  of sites after (if S_before is provided)

- S_before:

  Number of sites before, if you wish to keep this number constant in
  the analysis. Defaults to NULL, which makes S the same before and
  after.

- nB:

  Number of before measurements per site

- nA:

  Number of after measurements per site

- sd_pooled:

  Pooled standard deviation among all data points (on the transformed
  scale)

- baseline_mean:

  Mean of the response variable before the change in the original
  (untransformed) scale. Required for all typeTransform values to
  convert the detectable absolute delta back to a percent change.

- target_power:

  Desired power (default 0.8)

- alpha:

  Significance level (default 0.05)

- typeTransform:

  Character indicating the transformation applied to the response
  variable before analysis. One of `"none"`, `"log"`, `"sqrt"`, or
  `"asin"` for asin(sqrt) (default `"none"`).

- addValue:

  Value added to the response variable before transforming, to avoid
  issues with zeros (default 0). Ignored for `"asin"`.

## Value

The minimum detectable percentage change (e.g., 30 for a 30% change)
