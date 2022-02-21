# convenience functions


#' Robust natural log
#'
#' @param x Numeric vector or scalar
#'
#' @return natural log transformation of `x`, `Ln(0) ==0`
#' @export
#'
#' @examples
#' Ln(5)== log(5)
#' Ln(0)
#' # base log(0) returns unhelpful -Inf
#' log(0)
#'
Ln = function(x) {ifelse(x != 0, log(x), 0)}

#' Exponentiate `x` to power `y`
#'
#' @param x Numeric vector or scalar
#' @param y Numeric scalar
#'
#' @return
#' @export
#'
#' @examples
#'
#' all.equal(Exp(exp(1), 5),  exp(5) )
#' Ln(Exp(exp(1), 10)) ==10
Exp = function(x, y) {ifelse(x != 0, x^y, 0)}


#' Root mean squared log error
#'
#' @param x positive numeric vecotr
#' @param true_x numeric scalar, target x value, `mean(Ln(x))` by default
#'
#' @return scalar, root mean of the squared log ratio of observations to target value
#' @export
#'
#' @examples
#'
#' x <- rlnorm(n = 1e5)
#' rmsle(x) # should be very close to sdlog of lnorm, 1
#' y <- rnorm(n = 5000, mean = 12 )
#' rmsle(y, true_x = Ln(12))
#' 1/12

rmsle <- function(x, true_x = mean(Ln(x))){
  sqrt(mean((Ln(x) - true_x)^2))
}


#' Not in
#' @name %ni%
#' @export
#' @usage lhs \%ni\% rhs
#' @examples
#'
#' 6 %ni% 1:5

`%ni%` <- Negate(`%in%`)


