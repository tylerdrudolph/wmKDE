ReadMe
================
Tyler D. Rudolph

## **wmKDE: Weighted Mean Kernel Density Estimation**

### Description

Estimation of the weighted and un-weighted bivariate kernel probability
density function for spatially explicit[^1] population-level processes
given repeated measures of individual sample units. May be used to
estimate the utilization distribution of a wildlife or fish population
monitored via satellite (e.g. GPS collar transmitters) or acoustic
(e.g. PTT tags & receptor stations) data telemetry while accounting for
individual variation in point sample size, dispersion and quality,
thereby generating unbiased estimates of population distribution.

Two weighting strata are possible:

1.  Individual sample observations (e.g. GPS relocations) may be
    weighted to reflect variation in positional accuracy (e.g. PDOP,
    HDOP) and/or datum quality;

2.  Individual kernel distributions may be weighted to reflect variation
    in input data sufficiency, representativity and/or reliability.

## Installation

``` r
if(!require('devtools')) {
  install.packages('devtools')
}

devtools::install_github('https://github.com/tylerdrudolph/wmKDE')
```

## Demonstration

[^1]: Although current implementation requires a geo-spatially
    referenced vector feature input (sf object), this method can (and
    may eventually) be generalized to any comparable bivariate
    distribution problem.
