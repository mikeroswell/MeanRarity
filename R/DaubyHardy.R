#' Worker to iteratively solve ENS transformation of Hurlbert diversity
#'
#' @param s_k numeric Hurlbert rarefaction
#' @param ens_k_old numeric ENS estimate from previous iteration
#' @param k integer sample size parameter for rarefaction
#'
#' @return numeric scalar, updated ENS estimate
#'
ENSer <- function(s_k, ens_k_old, k){
  e_k_i = s_k/(1-(1-(1/ens_k_old))^k)
}

#' Hurlbert rarefaction
#'
#' Provides Hurlberts (1971) rarefied richness.
#'
#' @seealso vegan::rarefy()
#'
#' @template ab_template
#' @param k integer sample size parameter for rarefaction
#'
#' @return numeric scalar, unbiased estimate of rarefied richness, assuming
#'   abundances are sampled without replacement from a finite species pool
#'
#' @references
#' \insertRef{Hurlbert1971}{MeanRarity}
#'
#' @concept Estimators
#' @export
#'
#' @examples
#' hRare(1:10, 2)

hRare <- function(ab, k){
  summand = unlist(sapply(ab, function(x){
    1-choose(sum(ab)-x, k)/choose(sum(ab), k)
  }))
  return(sum(summand))
}


#' Estimate ENS based on Hurlbert rarefaction
#'
#' Implements estimate described in Dauby and Hardy 2011 for a class of
#' rarefaction-based ENS diversity estimates. These estimates suffer from
#' minimal bias and are quite efficient, while retaining some of the nice
#' properties of Hill diveristy metrics. They are parameterized by sample size
#' `k`, and when `k == 2` they are equivalent to Hill-Simpson diversity. One
#' interpretation is that this ENS is the number of species in a perfectly even
#' assemblage that would have the same rarefied richness as the focal
#' assemblage/sample. Larger k values emphasize rare species, and as k
#' approaches community size the Hulbert ENS approaches true richness. Unbiased
#' estimators are given for k < sample size.
#'
#' @template ab_template
#' @param k integer sample size parameter for rarefaction
#' @param maxit integer, maximum number of iterations
#' @param tol numeric, threshhold for convergence
#'
#' @concept Estimators
#' @return Numeric scalar: estimated Hurlbert ENS
#'
#' @references
#' \insertRef{Dauby2012}{MeanRarity}
#' \insertRef{Hurlbert1971}{MeanRarity}
#' @export
#'
#' @examples
#' ab = sample(10:50, 50, replace =TRUE)
#' hurl(ab, 2)
#' # not run
#' # hurl(ab, 1e5) # returns an error
#'
hurl <- function(ab, k, maxit = 1e5, tol = 1e-12){
  if(!all(ab == floor(ab))){stop("abundance vector must be integer values")}
  if(k < 1){stop("undefined for samples <1 individual")}
  if(k>sum(ab)){stop("k must be smaller than total abundance")}
  # if(sum(ab)<100){warning("with small abundances, sample size correction is severe")}
  i = 1
  s_k = hRare(ab, k) # for i = 0
  e_k_i = ENSer(s_k, s_k, k = k) # gives for i = 1
  while(i < maxit){
    i = i+1
    e_temp = e_k_i
    e_k_i = ENSer(s_k, e_temp, k =k)
    if(abs(e_k_i - e_temp) < tol) { break}
  }
  return(e_k_i)
}
