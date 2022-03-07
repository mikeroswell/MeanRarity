#' Ratio labels
#'
#'
#' @param x Numeric, on a logarithmic scale
#' @param base Scalar, base of the logarithm used
#'
#' @return String to be evaluated with \code{str2expression}
#' @export
#'
#' @examples
#' print_operator(((-1:3)))
#'
print_operator <- function(x, base = exp(1)){
  ifelse(sign(x) == -1,
         paste("NULL","%/%", base^abs(x))
         , ifelse(sign(x) == 1,
                  paste("NULL", "%*%", base^abs(x))
                  , paste0("bold(", base^x, ")")
                  )
  )
}


#' Compute breaks for ratio scale
#'
#' Function to compute tick marks evenly spaced on the log scale but with pretty
#' numbers on the ratio scale
#'
#' @param base Scalar, base of the logarithm in use (not implemented)
#' @param n scalar, target number of breaks (not implemented)
#'
#' @return Function to apply over a vector values to generate axis breaks
#' @export
#'
#' @examples
#' ggplot2::ggplot(data = data.frame(x= 1:6, y = seq(-1, 1.5, 0.5))
#'        , ggplot2::aes(x, y))+
#'              ggplot2::geom_point()+
#'              ggplot2::scale_y_continuous(
#'              , breaks = rat_breaks()
#'              ) +
#'              ggplot2::geom_hline(yintercept = 0, size = 0.2)
rat_breaks <- function( base = exp(1), n = 5){
  function(x){
    largest_integer = floor(exp(max(abs(x))))
    if(largest_integer >=2){
      one_side = Ln(floor(exp(-max(abs(x)): max(abs(x)))))
      both_sides = unique(c(one_side, -one_side))
      trun = c(-log(2)
               , both_sides[ifelse(min(x)<0, min(x), -1) < both_sides &
                              both_sides < ifelse(max(both_sides) >0
                                                  , 1.1 *max(both_sides), 2)]
               , log(2))
      br = trun}
    if(largest_integer <2){
      z = 0
      while(floor(10^z*(exp(max(abs(x)))-1))<2){
        z = z+1

      }
      if(-sign(max(x))==sign(min(x))){
        br = Ln(c( 1, 1+15*10^-z, 1+5*10^-z, 1/ (1+15*10^-z), 1/(1+5*10^-z)))
      }
      else{br = sign(max(x)) *
        Ln(c( 1, 1+5*10^-z, 1+2*10^-z, 1+10^-z, 1+5*10^-(z+1)))}
      # not sure it will see these little numbers
    }
    return(br)
  }
}


# # y-values are already on the log scale
# theme_set(theme_classic())
#
# #big healthy range
# ggplot(data = data.frame(x= 1:6, y = seq(-1, 4, 1))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#      labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
# # more reasonable values e.g. for our purposes
# ggplot(data = data.frame(x= 1:6, y = seq(-1, 1.5, 0.5))
#        , aes(x, y))+
#          geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
# # no integers
# ggplot(data = data.frame(x= 1:6, y = seq(-0.1, 0.15, 0.05))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
# # tiny changes
# ggplot(data = data.frame(x= 1:6, y = seq(-0.01, 0.015, 0.005))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
# # stay positive
# ggplot(data = data.frame(x= 1:6, y = seq(0.5, 3, 0.5))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
# # stay negative
# ggplot(data = data.frame(x= 1:6, y = seq(-3, -0.5, 0.5))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
#
#
# # big numbers cause issue, require stupid number crunching
# ggplot(data = data.frame(x= 1:6, y = seq(1, 11, 2))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)
#
# # not good for small numbers on one side of 0
# ggplot(data = data.frame(x= 1:6, y = seq(0.00001, 0.00006, 0.00001))
#        , aes(x, y))+
#   geom_point()+
#   scale_y_continuous(
#     labels = function(x){str2expression(print_operator(x))}
#     , breaks = rat_breaks()
#   ) +
#   geom_hline(yintercept = 0, size = 0.2)



