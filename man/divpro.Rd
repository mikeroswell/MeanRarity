% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DFuns.R
\name{divpro}
\alias{divpro}
\title{Generate Hill diversity profile}
\usage{
divpro(ab, ell_low = -1, ell_hi = 1, by = 0.001)
}
\arguments{
\item{ab}{A numeric vector of species abundances or relative abundances.}

\item{ell_low}{Scalar, minimum scaling exponent for diversity profile.}

\item{ell_hi}{Scalar, maximum scaling exponent for diversity profile.}

\item{by}{Scalar, size of step along scaling exponent continuum.}
}
\value{
Dataframe with the scaling exponent \{code{ell} and corresponding
Hill diversity \code{d}
}
\description{
Compute observed Hill diversity profile based on an abundance vector over a
range of scaling exponent values.
}
\details{
Hill diversity can be viewed as a continuous function of the scaling exponent
\eqn{\ell}{"ell"} and the relative abundance distribution. As \eqn{\ell}{"ell} increases,
so does the emphasis on rare species. It is traditional to view the profile
across \eqn{\ell = [-1, 1]}{"ell" = [-1, 1]} or \eqn{\ell = [-2, 1]}{"ell" = [-2, 1]}, and other authors
have visualized this with low values of \eqn{\ell}{"ell"} at the right instead of
left.
}
\examples{
divpro(c(20,8,5,4,2,1))
}
