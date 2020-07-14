# fucntions to Generate a semi-parametric SAD based on richness, and Hill-Simpson.


# define variables but also give them values for messing with
# totAb<-1e7 #total true abundance in community
# rich<-50 # true richness
# simpson<-40 #true simpson


#define my gamma distribution

#' Relative abundances given gamma distribution
#'
#' Wrapper for \code{\link[stats]{qgamma}} that
#' takes the number of species and a shape parameter
#' and gives relative abundance estimates for each species.
#'
#' @param x Shape parameter for a gamma distribution, a scalar
#' @param rich Total number of species in the SAD, an integer
#'
#' @seealso \code{\link[stats]{qgamma}, \link{fit_SAD}}
#'
#' @export
#' @examples
#' divers_gamma(rich = 10, x = 0.2)
#'
#' divers_gamma(rich = 10, x = 2)
divers_gamma<-function(rich, x){
    stats::qgamma(seq((1/rich)/2, 1-(1/rich)/2, (1/rich)), shape=x, scale=10/x)
}


#' Difference between target and realized diversity of simulated SAD (gamma
#' distribution).
#'
#' Subtracts realized inverse simpson diversity of the simulated species
#' abundance distribution from the target value. When this difference = 0 the
#' shape parameter of the gamma distribution is considered optimal.
#'
#' @param x Shape parameter for a either lognormal or gamma distribution, a scalar.
#' @param rich Total number of species in the SAD, an integer.
#' @param simpson Target value for inverse Simpson's concentration of the simulated
#'   SAD, a scalar.
#' @param distr Distribution type (currently "\code{lnorm"} or \code{"gamma"}) to call for
#'   "divers_" function, a character string.
#' @param totAB Not implemented, could have a finite-size version with a fixed #
#'   of individuals in pool.
#'
#' @return a scalar, the difference between empirical and target Hill-Simpson
#'   diversity
#'
#' @seealso \code{\link{dfun}}
#'
#' @export
#' @examples
#'
#' ur_distr(x = 2, rich = 10, simpson = 6, distr = "gamma")
#' ur_distr(x = 0.2, rich = 10, simpson = 6, distr = "gamma")
#' ur_distr(x = 1.235383382 , rich = 10
#'      , simpson = 6, distr = "gamma") #very close, depends on sensitivity.
ur_distr<-function(x,rich=rich, simpson=simpson, distr=distr, totAb=totAb, ...){
    simpson-dfun(get(paste0("divers_", distr))(rich, x), -1)}


#' Relative abundances given lognormal distribution
#'
#' This is a wrapper for \code{qlnrom} that
#'      takes the number of species and a shape parameter (sd of the log(lognormal)))
#'      and gives relative abundance estimates for each species.
#' @param x Shape parameter \code{sdlog}, a scalar.
#' @param rich Total number of species in the SAD, an integer.
#'
#' @return Numeric vector of species relative abundances
#' @seealso \code{\link[stats]{qlnorm}}
#'
#' @export
#' @examples
#' divers_lnorm(rich = 10, x = 1)
divers_lnorm<-function(rich, x){
    stats::qlnorm(seq((1/rich)/2, 1-(1/rich)/2, (1/rich)), meanlog=10, sdlog=x)
    }



#' Fit a SAD given diversity and a distributional assumption
#'
#' Takes a true richness, an inverse-Simpson diversity, and a distributional
#' assumption. Fits an optimal SAD given these constraints using
#' \code{stats::uniroot}.
#'
#' The way the SAD is fit, we assume infinite abundance, but finite and
#'      fixed diversity. In particular, richness is fixed. The parametric
#'      fits use continuous distributions. Species abundances are given at
#'      evenly spaced intervals along those continuous distributions. The distributions
#'      are described by a shape parameter (the scale parameter of the lognormal
#'      is fixed), and the value of that parameter is chosen such that the relative
#'      abundances assigned to the species give the target Hill-Simpson diversity.
#'
#'
#' @param rich Total number of species in the SAD, an integer.
#' @param simpson Hill-Simpson diversity of the SAD, a real number in [1,rich].
#' @param int_lwr Lower bound of search space for uniroot; default is a small
#'   number close to 0 to deal with potential boundary issues for
#'   \code{stats::uniroot} (a scalar).
#' @param int_uppr Upper bound of search space for \code{uniroot} (scalar).
#' @param distr Name of the distribution (\code{"lnorm"} or \code{"gamma"}).
#' @param totAb Not implemented; would be a necessary constraint for fitting
#'   finite communities.
#'
#' @return A list with three elements.
#'
#'   \code{distribution_info} Contains the name of the distribution and the fitted
#'   shape parameter value.
#'
#'   \code{community_info} gives richness, Hill-Shannon, and Hill=Simpson
#'   diversity of the SAD.
#'
#'   \code{rel_abundances} is a vector of relative abundances for each species
#'   in SAD.
#'
#'
#'
#' @seealso \code{\link[stats]{uniroot}}, \code{\link{dfun}}
#' @export
#' @examples
#' #works
#' fit_SAD(distr = "gamma") #works
#' fit_SAD(distr = "lnorm", rich = 50, simpson = 2) #works
#' fit_SAD(distr = "gamma", rich = 50, simpson = 2) #works
#' fit_SAD(distr = "gamma", rich = 10, simpson = 6)
#'
#' \dontrun{
#' fit_SAD(distr = "nonsense")  #returns error
#' fit_SAD(distr = "lnorm", rich = 50, simpson = 90) #returns error
#' }

fit_SAD<-function(rich = 50
                  , simpson = 40
                  , distr = "lnorm"
                  , int_lwr = 1e-04
                  , int_uppr = 1e2
                  , totAb = 1e7){
    #check feasibility; these should be actual errors in future versions
    if(simpson>rich| simpson<1){
        stop("Hill-Simpson diversity cannot be greater than richness nor less than 1")}

    #check dstr makes sense
    ifelse( !(distr %in% c("lnorm", "gamma")),
            stop("distr must be either `lnorm` or `gamma`"),

        #generate SAD optimized with uniroot to find x when ur_distr==0
                {fit_par = tryCatch(stats::uniroot(function(x){ur_distr(simpson=simpson, rich=rich, x=x, distr=distr, totAb=totAb)}
                                                , lower=int_lwr, upper=int_uppr)
                                 , error=function(e) message("test int_lwr and int_uppr as values for `x` in ur_distr(); output must have opposite signs")
                                 )

                #make sure to return rel abundances!
                abus=tryCatch(divers_lnorm(rich, fit_par$root)
                              , error=function(e) {message("did not fit param")
                    return(rep(100, length(rich)))
                                  }
                )

        # relative abundances
                abus<-abus/sum(abus)
        #return Hill-Shannon also
        shannon=dfun(abus, 0)
        # if(sum(abus==0)>0) print("WARNING: you simulated species with 0 abundance")

        return(list("distribution_info"=c("distribution"=distr, "fitted parameter"=fit_par$root)
                    , "community_info"=c("richness"=rich, "Hill-Shannon"=shannon, "Hill-Simpson"=simpson)
                    , "rel_abundances"=abus

                   ))}
    )
}





