# find_min_sites_analytical

Function to use power.t.test to calculate power for a paired t-test on
the means across

Function to use power.t.test to calculate power for a paired t-test on
the means across

## Usage

``` r
find_min_sites_analytical(
  S_grid,
  delta_target,
  sd_diff,
  distribution = c("normal", "nbinom", "binomial"),
  useTest = c("paired-t", "wilcoxon", "prop.test")
)

find_min_sites_analytical(
  S_grid,
  delta_target,
  sd_diff,
  distribution = c("normal", "nbinom", "binomial"),
  useTest = c("paired-t", "wilcoxon", "prop.test")
)
```

## Arguments

- delta_target:

  Target mean change

- sd_diff:

  Standard deviation of difference in means

- distribution:

  Distribution for simulated data. Must be "normal" here.

- useTest:

  Which test to use. Must be "paired-t" here.

- n_grid:

  Number of sites grid

## Value

Power for each number of sites in n_grid

Power for each number of sites in n_grid
