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
  logTransform = FALSE,
  logAdd = 0
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

- logTransform:

  Logical indicating whether to log-transform the response variable for
  estimating standard deviations on the log scale

- logAdd:

  Value to add to response variable before log-transforming to avoid
  issues with zeros (default 0)

## Value

A data frame containing the estimated within-site and between-site
standard deviations, the number of sites and before samples per site,
and mean proportion positive by site, summarized by group when groupVar
is provided.
