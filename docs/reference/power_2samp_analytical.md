# power_2samp_analytical

Analytical power for a two-sample (unpaired) t-test using the
non-central t distribution. Handles unbalanced group sizes.

## Usage

``` r
power_2samp_analytical(n_before, n_after, delta, sd_pooled, alpha)
```

## Arguments

- n_before:

  Number of before measurements (total across all sites)

- n_after:

  Number of after measurements (total across all sites)

- delta:

  Hypothesized mean change (absolute, on the scale of the response)

- alpha:

  Significance level (default 0.05)

- sd_w:

  Within-site standard deviation (or pooled standard deviation across
  all samples)
