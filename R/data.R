#' Bee specimens from 4 sites sampled with equal effort
#'
#' A dataset containing the bee species identities of all specimens
#' collected at each site during 15 30-minute transects.
#'
#' @format A tibble with 2019 rows and 7 variables:
#'  - *uniqueID*: unique specimen identifier
#'  - *bee*: concatenated Latin binomial for bee species ID
#'  - *start*: time of day that 30-minute sampling transect began, with
#'  arbitrary date value
#'  - *sr* concatenated site name and sampling round
#'  - *sday1* ordinal, day of sampling round (1, 2, or 3)
#'  - *sampling_round* integer between 1 and 5, indicating which 3-day sampling
#'  round at given site
"beeObs"

#' Bee specimens from 4 sites sampled with equal effort, summarized
#'
#' A dataset containing the number of specimens of each bee species
#' collected at each site during 15 30-minute transects.
#'
#' @format A tibble with 68 rows and 5 variables:
#'  - *bee*: concatenated Latin binomial for bee species ID
#'  - *\`Cold Soil_5\`*: number of records of given species from site "Cold
#'  Soil"  during sampling round 5
#'  - *\`Fox Hill_5\`*: number of records of given species from site "Fox Hill"
#'   during sampling round 5
#'  - *\`IAS_3\`*: number of records of given species from site "IAS" during
#'  sampling round 3
#'  - *\`Lord Stirling_4\`*: number of records of given species from site
#'  "Lord Stirling" during sampling round 5
"beeAbunds"

#' Rarefied diversity estimates for each site
#'
#' A dataset containing expected observed Hill diversity and estimated
#' asymptotic Hill diversity for each of the 4 community samples analyzed in
#' \insertCite{Roswell2020}{MeanRarity}.
#'
#' @references
#' \insertRef{Roswell2020}{MeanRarity}
#'
#' @format A tibble with 2015 rows and 8 variables:
#'  - *site*: Character, site ID
#'  - *individuals*: Integer, individuals in rarefied sample
#'  - *obsrich*: Numeric, mean observed species richness
#'  - *chaorich*: Numeric, mean estimated asymptotic richness from Chao1
#'  estimator
#'  - *obsshan*: Numeric, mean observed Hill-Shannon diversity
#'  - *chaoshan* Numeric, mean estimated asymptotic Hill-Shannon diversity
#'  - *obssimp* Numeric, mean observed Hill-Simpson diversity
#'  - *chaosimp* Numeric, mean estimated asymptotic Hill-Simpson diversity
"mean_ests"

#' Rarefied diversity estimates for each site, multiple standards
#'
#'?
#'
#' @seealso \insertRef{Chao2012a}{MeanRarity}.
#'
#' @format A tibble with 1008 rows and 6 variables:
#'  - *site*: Factor, site ID
#'  - *divind*: Factor, Hill diversity ("richness"
#'  , "Hill-Shannon", or "Hill-Simpson)
#'  - *etype* Factor, estimate type ("observed", or "asymptotic")
#'  - *diversity*: Numeric, mean estimated Hill diversity value
#'  - *method*: Factor, standardization method ("effort", "size", or "coverage")
#'  - *xax*: Numeric, number of transects in subsamples, mean number of
#'  individuals in subsamples, or mean sample coverage of subsamples
#'  - *obsshan*: Numeric, mean observed Hill-Shannon diversity
#'  - *chaoshan* Numeric, mean estimated asymptotic Hill-Shannon diversity
#'  - *obssimp* Numeric, mean observed Hill-Simpson diversity
#'  - *chaosimp* Numeric, mean estimated asymptotic Hill-Simpson diversity



"effort_rare"

#' Estimator accuracy and bias
#'
#' A simulated dataset with rmse and mean absolute error for different
#' approaches to estimating the true (asymptotic) Hill diversity based on finite
#' samples from simulated species abundance distributions. Hill diversity was
#' estimated with the naïve estimator, `rarity(sample_freq, l)`, the Chao and
#' Jost 2015 non-parametric asymptotic estimator `Chao_Hill_abu(sample_freq,
#' l)`, and God's estimator `GUE(sample_freq, true_freq, l)`.
#'
#' @references
#' \insertRef{Chao2015}{MeanRarity}
#'
#' @format A data.frame with 2376 rows and 8 variables:


#'  - *ell*: Scalar, scaling parameter for Hill diversity
#'  - *distribution* Factor, family, one of "lnorm" or "gamma"
#'  - *fitted.parameter* Factor, shape parameter for fitting SAD
#'  - *n*: Numeric, individuals sampled
#'  - *SAD*: index of SAD sampled
#'  - *estimator* Character, estimator whose error is described
#'  - *rmse* Numeric, root mean squared error of estimator
#'  - *bias* Numeric, mean difference between estimator and true diversity
#'

"errs"

