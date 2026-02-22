# power_for_percent_change

Function to calculate power for a given percent change in the mean
before vs. after change using power.t.test (paired t-test on site-level
differences).

## Usage

``` r
power_for_percent_change(
  percent_change,
  baseline_mean = NULL,
  nA,
  S,
  nB,
  sd_within,
  sd_delta,
  alpha,
  logTransform = FALSE,
  logAdd = 0
)
```

## Arguments

- percent_change:

  Target percent change in the original response variable

- baseline_mean:

  Mean of the original response variable before the change

- nA:

  Number of after samples per site

- S:

  Number of sites

- nB:

  Number of before samples per site

- sd_within:

  Standard deviation within sites (log calculated if logTransform=TRUE)

- sd_delta:

  Standard deviation of true changes among sites

- alpha:

  Significance level

- logTransform:

  Logical indicating whether to work on log(y + logAdd) scale (default
  FALSE)

- logAdd:

  Value to add to response variable before log-transforming to avoid
  issues

## Value

Power for the given percent change
