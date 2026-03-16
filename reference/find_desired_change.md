# find_desired_change

Function to convert a desired change in the original response variable,
specified as either an absolute change or a percent change, to the
corresponding change on the scale of the response variable (e.g. log or
sqrt) if a transformation is specified, using the baseline mean to
convert percent changes to absolute changes on the original scale before
applying the transformation.

## Usage

``` r
find_desired_change(
  change_type,
  change_value,
  baseline_mean,
  typeTransform = c("none", "log", "sqrt"),
  addValue = 0
)
```

## Arguments

- change_type:

  Character indicating whether the desired change is specified as an
  "absolute" change or a "percent" change in the original response
  variable.

- change_value:

  Numeric value of the desired change in the original response. For
  example, if change_type is "percent", a value of 30 means a desired
  change of 30% (i.e. multiply original mean by 1.3).

- baseline_mean:

  Mean of the original response variable before the change

- typeTransform:

  Character indicating the transformation to apply to the response
  variable, one of "none", "log", or "sqrt" (default "none").

- addValue:

  Value to add to response variable before transforming (default 0).
