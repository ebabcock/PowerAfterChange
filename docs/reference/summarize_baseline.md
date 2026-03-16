# summarize_baseline

Summarize baseline data to estimate within-site and between-site
standard deviations, and optionally log-transform the response variable
to estimate these on the log scale.

## Usage

``` r
summarize_baseline(
  baseline,
  siteVar = "site",
  responseVar = "y",
  groupVar = NULL,
  typeTransform = "none",
  addValue = 0
)
```

## Arguments

- baseline:

  Data frame containing baseline measurements

- siteVar:

  Name of the site variable in baseline data

- responseVar:

  Name of the response variable in baseline data

- groupVar:

  Optional character vector of grouping variables to summarize by

- addValue:

  Value to add to response variable before transforming (default 0)

- typeTansform:

  Character indicating transformation the response variable,
  "none","log", or "sqrt" (default "none")"

## Value

A data frame containing the estimated mean and standard deviation of the
response variable across all sites and data points
(grand_mean,grand_sd), and the within and across site standard
deviations (sd_within, sd_between) in the original scale, and optionally
the log or sqrt transformed scale if typeTransform is "log" or "sqrt",
with an added constant if provide d(addValue). The number of sites in
the baseline is provided (nB), or, if the number of samples per site is
not consistent, the number of complete and the number of sites with
fewer samples (nb_complete and nB_incomplete). The proportion of
positive values in the response variable is also calculated
(prop_positive) to help assess whether a log transformation is
appropriate. If groupVar is provided, these statistics are calculated
separately for each group defined by groupVar.
