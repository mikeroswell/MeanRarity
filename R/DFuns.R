# Here are the key mean rarity functions

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
pfun=function(x, pow){
  if (pow==0) return(log(x))
  r <- sign(pow)*(x)^pow
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
ipfun=function(x, pow){
  if (pow==0) return(exp(x))
  x<-ifelse(sign(pow)*x<0,0,x) #added so that ggplot padding doesnt introduce negative values to scale
  r <- (sign(pow)*(x))^(1/pow)
  return(r)
}


#' Estimate Hill diversity: the mean species rarity
#'
#' We parameterize Hill diversity \eqn{D} as a the frequency-weighted mean species rarity, with scaling exponent l
#'      \deqn{D = \sum{p_i * r_i^{\ell}}^{-\ell}}
#'      where rarity of species i \eqn{r_1 = 1/p_i}.
#'      When \eqn{\ell = 0} this is defined base on the limit from the left and the right, which is the
#'      geometric mean \deqn{\exp(\frac{\sum{p_i * \ln(r_i)}} \sum{p_i}})
#'
#' @details This is equivalent to the \eqn{q} notation of Jost 2006
#'      \deqn{D=\sum{p_i^q}^{1-q}}
#'      where \eqn{q=1-l}.
#'
#'
#' @param ab A numeric vector of species abundances or relative abundances.
#' @param l Scaling exponent for the mean, can be any real number.
#'
#' @return Generalized mean community rarity with scaling exponent \code{"l"}.
#'
#' When \code{l = 1}, arithmetic mean rarity (species richness).
#'
#' When \code{l = 0}, geometric mean rarity (Hill-Shannon diversity), Shannon's entropy
#'    (Shannon and Weaver 1963) exponentiated.
#'
#' When \code{l = -1}, harmonic mean rarity (Hill-Simpson diversity,
#'    the inverse of the Simpson concentration (Simpson 1949))
#'
#' @seealso \code{\link{pfun}}, \code{\link{ipfun}}
#'
#' @export
#' @examples dfun(c(20,8,5,4,2,1), -1)
dfun<-function(ab, l){
  ab<-ab[ab!=0]
  rp <- ab/sum(ab)
  if(l==0) {return(exp(sum(rp*log(1/rp))))}
  return(ipfun(sum(rp*pfun(1/rp, l)),l)) #removed potentially problematic sign corrections
  # return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/rp, l)),l)) # is it possible there was a mistake here?
}
