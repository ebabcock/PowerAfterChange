# find_n_after_2samp

Find the minimum number of after measurements per site needed to achieve
target power for a two-sample (unpaired) t-test that pools all before
and after measurements, ignoring within-site correlation.

## Usage

``` r
find_n_after_2samp(
  S,
  S_before = NULL,
  nB,
  delta,
  sd_pooled,
  target_power = 0.8,
  alpha = 0.05,
  n_grid = 1:50,
  typeTransform = c("none", "log", "sqrt", "arcsin"),
  addValue = 0,
  baseline_mean = NULL
)
```

## Arguments

- S:

  Number of sites

- S_before:

  Number of sites before, if you wish to keep this number constant in
  the analysis. Defaults to NULL, which keeps S sites before and after.

- nB:

  Number of before measurements per site

- delta:

  Hypothesized mean change (absolute, on the original untransformed
  scale of the response). When `typeTransform` is not `"none"`, this
  value is converted to the transformed scale internally using
  `find_desired_change` before computing power.

- sd_pooled:

  Standard deviation among all samples (on the transformed scale)

- target_power:

  Desired power level (default 0.8)

- alpha:

  Significance level (default 0.05)

- n_grid:

  Grid of after measurements to evaluate (default 1:50)

- typeTransform:

  Character indicating the transformation applied to the response
  variable before analysis. One of `"none"`, `"log"`, `"sqrt"`, or
  `"arcsin"` for arcsin(sqrt) (default `"none"`).

- addValue:

  Value added to the response variable before transforming, to avoid
  issues with zeros (default 0). Ignored for `"arcsin"`.

- baseline_mean:

  Mean of the original (untransformed) response variable before the
  change. Required when `typeTransform` is not `"none"` in order to
  convert `delta` to the transformed scale.

## Value

A list with `n_star` (minimum after measurements for target power) and
`curve` (data frame of n_after and power)
