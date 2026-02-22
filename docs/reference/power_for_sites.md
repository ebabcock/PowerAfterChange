# power_for_sites

Function to estimate power for a paired t-test on the average of a given
number of sites before and after a changepoint.

## Usage

``` r
power_for_sites(
  S,
  nB,
  nA,
  delta,
  sd_w,
  sd_d = 0,
  alpha = 0.05,
  nsim = 2000,
  seed = 1,
  distribution = c("normal", "nbinom", "binomial"),
  useTest = c("paired-t", "wilcoxon", "prop.test"),
  nbinom_mu = NULL,
  nbinom_disp = NULL,
  binomial_size = NULL,
  binomial_prob = NULL
)
```

## Arguments

- S:

  Number of sites

- nB:

  Number of before measurements per site

- nA:

  Number of after measurements per site

- delta:

  Hypothesized mean change

- sd_w:

  Within-site standard deviation

- sd_d:

  Between-site standard deviation of true changes

- alpha:

  Significance level, defaults to 0.05

- nsim:

  Number of simulations to run

- seed:

  Random seed for reproducibility

- distribution:

  Distribution for simulated data. One of "normal", "nbinom", or
  "binomial".

- useTest:

  Which test to use. One of "paired-t", "wilcoxon", or "prop.test".

- nbinom_mu:

  Mean parameter for negative binomial (mu)

- nbinom_disp:

  Dispersion (size) parameter for negative binomial

- binomial_size:

  Size parameter (trials) for binomial

- binomial_prob:

  Probability parameter for binomial

## Value

Estimated power for the given number of sites
