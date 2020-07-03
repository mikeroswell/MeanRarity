###### Here are the key mean rarity functions

###############################################################
##############################################################

#transformation and back-transformation functions

#' power transformation
#' @export
pfun=function(x, pow){
  if (pow==0) return(log(x))
  r <- sign(pow)*(x)^pow
  return(r)
}

#' inverse function for power transformation
#' @export
ipfun=function(x, pow){
  if (pow==0) return(exp(x))
  x<-ifelse(sign(pow)*x<0,0,x) #added so that ggplot padding doesnt introduce negative values to scale
  r <- (sign(pow)*(x))^(1/pow)

  return(r)
}

#' estimate Hill diversity
#'
#' We parameterize Hill diversity \eqn{D} as a the frequency-weighted mean species rarity, with scaling exponent l
#' \deqn{D=\sum{p_i*r_i^{l}}^{-l}}
#' where rarity of species i \eqn{r_1=1/p_i}
#'
#' This is equivalent to the "q" notation of Jost 2006
#' \deqn{D=\sum{p_i{q}}^{1-q}}
#' where \eqn{q=1-l}
#'
#' @param ab a numeric vector of species abundances or relative abundances
#' @param l scaling exponent for the mean, can be any real number
#'
#' @return generalized mean community rarity with scaling exponent "l".
#'
#' When \code{l = 1}, arithmetic mean rarity (species richness).
#' When \code{l = 0], geometric mean rarity (Hill-Shannon diversity)
#' When \code{l = -1}, harmonic mean rarity (Hill-Simpson diversity, the inverse of the Simpson concentration (Simpson 1949))
#'

#' @export
dfun<-function(ab, l){
  ab<-ab[ab!=0]
  rp <- ab/sum(ab)
  if(l==0) {return(exp(sum(rp*log(1/rp))))}
  return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/rp, l)),l))
}
