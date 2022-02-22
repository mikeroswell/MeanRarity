# MeanRarity
You can install the R package MeanRarity with the following command:

`devtools::install_github("mikeroswell/MeanRarity")`

This R package allows estimating and visualizing Hill diversity (aka "Hill
numbers") in terms of the average species "rarity."

The average species rarity may be simpler to interpret than a transformed,
generalized entropy. The rarity parameterization clarifies the role of the
scaling exponent, and this R package provides illustration tools to visualize
that scaling concretely. Furthermore, this framework may provide an entry point
for developing and testing diversity estimators and their related interval
estimates (e.g. confidence intervals); this package contains novel tools to
simulate species abundance distributions that should streamline that process.

## Three key functions
`rarity()` computes Hill diversity of a sample, parameterized in terms of the
mean rarity. 

`rarity_plot()` makes balance plots to visualize mean rarity.

`fit_SAD()` simulates a species abundance distribution with known richness,
Hill-Simpson diversity, and a distributional assumption (lognormal or gamma).

## Details
`rarity(x, l)` is equivalent to `vegan::renyi(x, scale = 1-l, hill = T)`. 
Re-parameterizing Hill diversity in terms of species rarity can provide 
conceptual clarity.

"Balance plots" called with `rarity_plot()` illustrate the computation in 
`rarity()`. 

`fit_SAD()` is a potentially novel way to simulate abundance distributions that
always gives precisely the same relative abundances for a given set of
parameters, and may be a helpful tool for exploring the behavior of diversity
metrics under different sampling scenarios.

## Vignettes
We are working on vignettes now. Currently Estimating_Mean_Rarity provides code
to make the figures in Roswell et al. 2021 Oikos
https://doi.org/10.1111/oik.07202, Gods_estimator explores a playful extension
of the tools in this package, and Using_MeanRarity is the beginning of a more
traditional "user's guide" type vignette. We have some work to do on this; feel
free to share suggestions and questions!

Thank you!

<!-- badges: start -->
  [![R-CMD-check](https://github.com/mikeroswell/MeanRarity/workflows/R-CMD-check/badge.svg)](https://github.com/mikeroswell/MeanRarity/actions)
<!-- badges: end -->
