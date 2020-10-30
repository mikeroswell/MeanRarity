# MeanRarity
Repository for development of MeanRarity R package

The goals of the package are to document and share the code and, perhaps more importantly, quantitative "methods" for Roswell et al. 2020 _Oikos_ and to provide a discrete, shareable codebase for working with mean rarity and associated CI.

In the Author's estimation, the cool functions are `rarity()`, which simply
computes the Hill diversity, but parameterizes this in terms of the mean rarity,
`rarity_plot()`, which makes balance plots to visualize mean rarity, and
`fit_SAD()`, which is a potentially novel way to simulate species abundance
distributions with known diversity, useful for testing diversity estimator
performance. Other documented functions support the use of these core tools or
are specific tools for analyses in the _Oikos_ MS.

TODO: 

 * look at the character that should be \eqn{\ell} use in vignette, causes errors
now.
 * Finish vignette(s)
   * will there be one for the MS and another more introducing the package, or
just one?
   * add/fix citations... see below, not comprehensive
 * Create tests
 * Add citation to Oikos paper
 * Add cites to Shannon and Weaver 1963, Simpson 1949, Chao and Jost 2015 MEE, 
Jost 2006, Patil and Taillie, Roswell et al 2019. 




