#funny estimators

#' Bad estimate of relative abundance with an incorrect finite-size correction
#'
#' @param x A numeric vector of species abundances
#'
#' @return An incorrect estimate of relative abundance
#'
#' @export
fsr<-function(x){
  (x-0.8)/(sum(x)-0.8)
}


#' Incorrect diversity estimator based on finite-size fudge
#'
#' This is a riff on \code{\link{dfun}} that uses \code{\link{fsr}}, which is not the correct way to deal with finite
#' size corrections, though it is similar to E.H. Simpson's estimator. Rather than subtracting 1, though, it subtracts 0.8
#'
#' @param ab Numeric vector of species abundances
#' @param l Exponent determining type of mean (scalar)
#'
#' @return Incorrect estimate of Hill Diversity with a goofy bias correction term (a scalar)
#'
#' @seealso \code{\link{dfun}}
#'
#' @export
fsd<-function(ab, l){
  ab<-ab[ab!=0]
  rp <- ab/sum(ab)
  fs<-fsr(ab)
  if(l==0) {return(exp(sum(rp*log(1/fs))))}
  return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/fs, l)),l))
}

#' Simpson's estimator of Hill-Simpson Diversity (we think?)
#'
#' E.H. Simpson's 1949 estimator of Simpson's Concentration is an unbiased estimator. It is not an unbiased Hill-Simpson estimator, but it is a pretty good one
#'
#' @param ab Numeric vector of species abundances
#'
#' @return the reciprocal of the Simpson's concentration, a low-bias estimator for Hill-Simpson diversity (scalar)
#' @export
sApp <- function(ab){
  n <- sum(ab)
  return(1/(
    sum(ab/n)*((ab-1)/(n-1)))
  ))
}
