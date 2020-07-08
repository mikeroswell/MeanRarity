#' MeanRarity: A package for computing Hill diversity under the "mean rarity"
#' framework.
#'
#' The MeanRarity package provides a suite of functions related to computing,
#' estimating, and visualizing Hill diversity: For references on mean rarity,
#' see Roswell et al. 2020 Oikos, and Patil and Taillie, and Jost 2006
#'
#' @section Computing Hill diversity: \code{\link{dfun}} computes Hill diversity
#'   as parameterized by Roswell et al. 2020, by computing the generalized
#'   weighted mean of species rarities (i.e. the reciprocal of relative
#'   abundance). The "link" functions for the generalized mean are
#'   \code{\link{pfun}} and its inverse \code{\link{ipfun}}
#'
#'
#' @section Visualizing Mean Rarity: The code for balance plots, also in the the
#'   web app meanrarity (URL) is also contained in this package
#'
#' @docType package
#' @name MeanRarity
