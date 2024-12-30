<!-- badges: start -->
[![R-CMD-check](https://github.com/olliemcdonald/bestdr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/olliemcdonald/bestdr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Bayesian Esitmation of Stochastic processes for Dose-Response (BESTDR)

BESTDR is an R package that uses Stan to estimate parameters for continuous time Markov
branching processes (CTMBP) approximated from a normal distribution. The package requires a user
to specify the model based on the transition probabilities of types along with a possible
"dose-response" statistical model connecting each rate to an underlying set of parameters.
The model is translated into Stan code that can be compiled and run.

This is a companion to the manuscript **Modeling mechanism-specific concentration response
in cell culture** available on
[https::/github.com/olliemcdonald/bestdr](https://github.com/olliemcdonald/bestdr).

### References:
McDonald, T.O., Bruno, S., Roney., J.P., Zervantonakis, I., Michor, F. BESTDR: Bayesian quantification of
mechanism-specific drug response in cell culture

## Installation:
Install devtools then install the package from Github.
``` r
install.packages('devtools')
devtools::install_github('olliemcdonald/bestdr')\
```

[CmdStan](https://github.com/stan-dev/cmdstan) needs to be installed as well. This can be done with the [cmdstanr](https://mc-stan.org/cmdstanr/) package.
``` r
cmdstanr::install_cmdstan()
```

## Approach
We developed methods for dose-response estimation of cell kinetic parameters along with scripts to create custom models as well as customize the dose-response parameter. We model cell growth as a multitype branching process and use the Central Limit Theorem for Branching Processes to approximate the likelihood, allowing us to easily estimate parameters by numerically solving the first two moments. These calculations are done using Stan and its internal ODE solvers. The Appendix contains an 'Extended Methods' section with more details abotu the mathematics.

## Getting Started
This package can be installed as an R package allowing easier customization of models. A vignette along with example script describe the process to create custom models. For basic concentration-response modeling of total cell counts, the birth and death rates can each be estimated. We've included a vignette for simulating and estimating the birth and death ratesfor a single cell type as well as multiple cell types in a hierarchical model that describes how data should be structured for estimation.

### For more details, look at the vignettes, manuscript, and accompanying supplementary information
