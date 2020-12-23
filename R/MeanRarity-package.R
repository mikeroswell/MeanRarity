#' MeanRarity: A package for computing Hill diversity under the "mean rarity"
#' framework.
#'
#' The MeanRarity package provides a suite of functions related to computing,
#' estimating, and visualizing Hill diversity: For references on mean rarity,
#' see \insertCite{Roswell2020}{MeanRarity}, \insertCite{Patil2012}{MeanRarity},
#' and \insertCite{Jost2006}{MeanRarity}
#'
#' @section Computing Hill diversity:
#'  \code{\link{rarity}} computes Hill
#'   diversity as parameterized by \insertCite{Roswell2020}{MeanRarity}, by
#'   computing the generalized weighted mean of species rarities (i.e. the
#'   reciprocal of relative abundance). The "link" functions for the generalized
#'   mean are \code{\link{pfun}} and its inverse \code{\link{ipfun}}
#'
#' @section Simulating Species Abundance Distributions:
#'   \insertCite{Roswell2020}{MeanRarity} developed a method to simulate a
#'   species abundance distribution (SAD) based on the \[Hill\] diversity of the
#'   assemblage. This method fits a discretization (based on quantile values) of
#'   a continuous parametric distribution (currently gamma or log-normal) based
#'   on the known, simulated richness and (currently) Hill-Simpson diversity.
#'   The function to do this is called \code{\link{fit_SAD}}
#'
#'
#' @section Visualizing Mean Rarity:
#' The code for rarity balance plots, also in
#'   the web app \url{https://mean-rarity.shinyapps.io/rshiny_app1/} is
#'   also contained in this package. \code{\link{rarity_plot}} makes rarity
#'   balance plots. \code{\link{radplot}} makes plots of rank-abundance
#'   distributions, with or without log-transforming abundances
#'
#' @importFrom rlang .data
#'
#' @references
#' \insertAllCited{}
#'
#' @keywords internal
#' @docType package
"_PACKAGE"
# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
