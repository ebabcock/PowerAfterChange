# power_for_nA

Function to calculate power for given number of after samples for a
paired t-test of the difference in means before vs. after change using
power.t.test

## Usage

``` r
power_for_nA_analytical(nA, S, nB, delta, sd_within, sd_delta, alpha)
```

## Arguments

- nA:

  Number of after samples per site

- S:

  Number of sites

- nB:

  Number of before samples per site

- delta:

  Hypothesized mean change

- sd_within:

  Standard deviation within sites

- sd_delta:

  Standard deviation of true changes among sites

- alpha:

  Significance level

## Value

Power for the given number of after sample
