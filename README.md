# MeanRarity
You can install the R package MeanRarity with the following command:

`devtools::install_github("https://github.com/mikeroswell/MeanRarity.git")`

This R package allows estimating and visualizing Hill diversity (aka hill numbers) in terms of the average species "rarity."

`rarity()` computes Hill diversity of a sample,  parameterized in terms of the mean rarity
`rarity_plot()` makes balance plots to visualize mean rarity
`fit_SAD()` simulates a species abundance distribution with known richness, Hill-Simpson diversity, and a distributional assumption (lognormal or gamma).

`rarity(x, l)` is equivalent to `vegan::renyi(x, scale = 1-l, hill = T)`. It is included to demonstrate how reparameterizing Hill diversity in terms of species rarity can provide conceptual clarity.

"Balance plots" called with `rarity_plot` illustrate the computation in `rarity()`. 

`fit_SAD()` is a potentially novel way to simulate abundance distributions that always gives precisely the same relative abundances for a given set of parameters, and may be a helpful tool for exploring the behavior of diversity metrics under different sampling scenarios.
