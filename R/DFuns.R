## Here are the key mean rarity functions

#' Power transformation
#'
#' Transform a vector by raising to a power.
#'
#' @param x A numeric vector.
#' @param pow The value of the exponent to which all values are raised.
#'
#' @return A numeric vector with same length as \code{x}.
#'
#' @export
pfun = function(x, pow){
  if (pow == 0) return(log(x))
  r <- sign(pow) * (x)^pow
  return(r)
}

#' Inverse function for power transformation
#'
#' Transform a vector by raising to a power... specifically 1/pow.
#'
#' @param x A numeric vector.
#' @param pow The reciprocal of the exponent to which all values
#' in \code{x} are raised.
#'
#' @return A numeric vector with same length as \code{x}.
#'
#' @export
ipfun = function(x, pow){
  if (pow == 0) return(exp(x))
  x <- ifelse(sign(pow) * x < 0, 0, x) #added so that ggplot padding doesn't introduce negative values to scale
  r <- (sign(pow) * (x))^(1 / pow)
  return(r)
}


#' Compute Hill diversity: the mean species rarity
#'
#' Compute the empirical Hill diversity from abundances or relative abundances.
#' Hill diversity is also the mean species rarity.
#'
#' We parameterize Hill diversity \eqn{D} as a the frequency-weighted mean
#' species rarity, with scaling exponent l \deqn{D = \sum{p_i *
#' r_i^{\ell}}^{-\ell}} where rarity of species i \eqn{r_1 = 1/p_i}. When
#' \eqn{\ell = 0} this is defined base on the limit from the left and the right,
#' which is the geometric mean \deqn{\exp(\frac{\sum{p_i * \ln(r_i)}}{\sum{p_i}})}
#'
#' This is equivalent to the \eqn{q} notation of Jost 2006
#'      \deqn{D=\sum{p_i^q}^{\frac{1}{1-q}}}
#'      where \eqn{q=1-l}.
#'
#' This function can also be called with \code{dfun()}
#'
#'
#' @template ab_template
#' @template l_template
#' @template q_template
#'
#' @return Generalized mean community rarity with scaling exponent \code{"l"}.
#'
#' When \code{l = 1}, arithmetic mean rarity (species richness).
#'
#' When \code{l = 0}, geometric mean rarity (Hill-Shannon diversity), Shannon's
#'    entropy \insertCite{Shannon1963}{MeanRarity} exponentiated.
#'
#' When \code{l = -1}, harmonic mean rarity (Hill-Simpson diversity),
#'    the inverse of the Simpson concentration
#'    \insertCite{Simpson1949}{MeanRarity}.
#'
#' @seealso \code{\link{pfun}}, \code{\link{ipfun}}
#' @references
#' \insertRef{Simpson1949}{MeanRarity}
#' \insertRef{Shannon1963}{MeanRarity}
#'
#' @concept Computation
#'
#' @export
#' @examples rarity(c(20,8,5,4,2,1), 1) #species richness
#' rarity(c(20,8,5,4,2,1), 0) # Hill-Shannon diversity
#' rarity(c(20,8,5,4,2,1), -1) # Hill-Simpson diversity
#' rarity(c(20,8,5,4,2,1), q = 2) # The parameter `q` can be used instead for
#' # traditional Hill number parameterization
rarity = function(ab, l, q = NULL){
  if(!is.null(q)){
    l = 1-q
    warning("l has been set to 1-q")
  }
  ab = ab[ab != 0]
  rp = ab/sum(ab)
  if(l == 0){return(exp(sum(rp * log(1/rp))))}
  D<-ipfun(sum(rp * pfun(1/rp, l)), l)
  if(!is.finite(D)){
    warning(
      "when `abs(l)` is very large, `rarity()` can return infinite results (multiplication errors with double-precision numbers). As l->Inf, rarity -> the maximum species rarity, ~N. As l -> -Inf, rarity -> the minimum species rarity, ~1."
    )
  }
  return(D) #removed potentially problematic sign corrections
  # return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/rp, l)),l))
  # is it possible there was a mistake here?
}

dfun <- rarity

#' Generate Hill diversity profile
#'
#' Compute observed Hill diversity profile based on an abundance vector over a
#' range of scaling exponent values.
#'
#' Hill diversity can be viewed as a continuous function of the scaling exponent
#' \eqn{\ell}{"ell"} and the relative abundance distribution. As \eqn{\ell}{"ell"} increases,
#' so does the emphasis on rare species. It is traditional to view the profile
#' across \eqn{\ell = [-1, 1]}{"ell" = [-1, 1]} or \eqn{\ell = [-2, 1]}{"ell" = [-2, 1]}, and other authors
#' have visualized this with low values of \eqn{\ell}{"ell"} at the right instead of
#' left.
#'
#' @template ab_template
#' @param ell_low Scalar, minimum scaling exponent for diversity profile.
#' @param ell_hi Scalar, maximum scaling exponent for diversity profile.
#' @param by Scalar, size of step along scaling exponent continuum.
#' @param use.q Logical, use traditional q parameterization, where q= 1-l.
#'
#' @return Dataframe with the scaling exponent \code{ell} and corresponding
#'   Hill diversity \code{d}
#'
#' @concept Computation
#'
#' @export
#' @examples divpro(c(20,8,5,4,2,1))
#' divpro(c(20,8,5,4,2,1), use.q =TRUE) # Option to use q = 1-ell for traditional
#' # Hill number parameterization

divpro <- function(ab, ell_low = -1, ell_hi = 1, by = 0.001, use.q = FALSE){
  ell = seq(ell_low, ell_hi, by = by)
  d = sapply(ell
         , function(l){rarity(ab, l)}
         )
  if(use.q){
    return(data.frame(q = 1-ell, d))
  }
  if(!use.q){
  return(data.frame(ell, d))
  }
}

#' Observed and asymptotic diversity
#'
#' Computes observed and asymptotic Hill diversity estimates from a vector of
#' integer abundances.
#'
#' Only integer abundances allowed in \code{ab}, as asymptotic estimator relies
#' on sampling theory for individuals.
#'
#' @return Dataframe with total abundance and observed and asymptotic diversity
#' estimates for the Pythagorean mean rarities: richness, exponentiated Shannon,
#' and inverse Simpson.
#'
#' @template ab_template
#'
#'
#' @seealso \code{\link{rarity}}, \code{\link{Chao_Hill_abu}}
#'
#' @export
#' @examples obs_est(c(20,8,5,4,2,1))

obs_est = function(ab){
  obsrich = rarity(ab, 1)
  obsshan = rarity(ab, 0)
  obssimp = rarity(ab, -1)
  chaorich = Chao_Hill_abu(ab, l = 1)
  chaoshan = Chao_Hill_abu(ab, l = 0)
  chaosimp = Chao_Hill_abu(ab, l = -1)
  coverage = Chat.Ind(ab)
  return(data.frame(
    n = sum(ab)
    , coverage
    , obsrich
    , chaorich
    , obsshan
    , chaoshan
    , obssimp
    , chaosimp
  ))
}

#' Normalized Hill evenness
#'
#' Computes assemblage evenness measure based on the ratio of Hill diversities.
#'
#' There are many valid ways to compute evenness, and this one is the
#' "normalized slope" of the Hill diversity profile. The normalization means
#' that for any richness, the maximum evenness is 1, and the minimum, achieved
#' when 1 species has practically all of the relative abundance in the
#' assemblage, is 0. This evenness metric was defined by Chao and Ricotta
#' (2019) as measure "E3". Like Hill diversity itself, evenness is considered a
#' function of scale, as well as the relative abundances in the assemblage.
#'
#' @seealso e2d
#' @concept Computation
#'
#' @template ab_template
#' @template l_template
#' @template q_template
#'
#' @export

e3Fun = function(ab, l, q = NULL){
  if(!is.null(q)){
    l = 1-q
    warning("l has been set to 1-q")
  }
  (rarity(ab, l) - 1)/(rarity(ab, 1)-1)
}

#' Hill diversity from evenness
#'
#' returns Hill diversity based on evenness and richness
#'
#' This function is the inverse of e3Fun
#'
#' @seealso e3Fun
#' @concept Computation
#'
#' @param e Scalar, evenness value between 0 and 1
#' @param rich Integer, species richness
#'
#' @export

e2d = function(e, rich){
  d = e*(rich-1)+1
  return(d)
}

