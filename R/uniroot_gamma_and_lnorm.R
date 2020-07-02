# Generates a semi-parametric SAD based on richness, and simpson, and total abundance. Requires a bit of thought b/c function requires a "prior" for the "shape" parameter of either "gamma" or "lnorm." For gamma, seems like if it tries the a shape parameter too close to 0, it returns garbage in the form of ur_gamma(0)==-Inf.

#returns a list. First element is SAD type and shape parameter. Second element is summary community info (total abundance, richness, shannon, simpson). Final element is abundance vector. 0 abundances throw warning but not error.


#define variables but also give them values for messing with
# totAb<-1e7 #total true abundance in community
# rich<-50 # true richness
# simpson<-40 #true simpson


#define my gamma distribution

#' relative abundances given gamma disribution
#'
#' This is a wrapper for `qgamma` that
#' takes the number of species and a shape parameter
#' and gives relative abundance estimates for each species
#' @param x Shape paramter for a gamma distribution, a scalar
#' @param rich Total number of species in the SAD, an integer
#'
#' @export
#' @examples
#' divers_gamma(rich = 10, x = 0.2)
#'
#' divers_gamma(rich = 10, x = 2)
divers_gamma<-function(rich, x){
    qgamma(seq((1/rich)/2, 1-(1/rich)/2, (1/rich)), shape=x, scale=10/x)
}



#function for uniroot to optimize, according to simpsons

#' difference between target and realized diversity of simulated SAD (gamma distribution)
#'
#' Subtracts realized inverse simpson diversity of
#' the simulated species abundance distribution from the target value.
#' When this difference = 0 the shape paramter of the gamma distribution is
#' considered optimal
#'
#' @param x Shape paramter for a gamma distribution, a scalar
#' @param rich Total number of species in the SAD, an integer
#' @param simpson Target value for inverse simpson diversity of the simulated SAD, a scalar
#' @param totAB Not implemented, could have a finite-size version with a fixed # of individuals in pool
#'
#' @export
#' @examples
#'
#' ur_gamma(x = 2, rich = 10, simpson = 6)
#' ur_gamma(x = 0.2, rich = 10, simpson = 6)
#' ur_gamma(x = 1.235383382 , rich = 10, simpson = 6) #very close, depends on sensitivity.
ur_gamma<-function(x,rich=rich, simpson=simpson, totAb=totAb,...){
    simpson-dfun(divers_gamma(rich, x), -1)}

# gamShape<-uniroot(ur_gamma, lower=1e-2, upper=1e2)
# gamShape

# ur_gamma(1e-5)

#define my lognormal distribution
divers_lnorm<-function(rich, x){
    qlnorm(seq((1/rich)/2, 1-(1/rich)/2, (1/rich)), meanlog=10, sdlog=x)
    }

#function for uniroot to optimize, according to simpsons
ur_lognorm<-function(x, rich=rich, simpson=simpson, ...){ ddiff=simpson-dfun(divers_lnorm(rich, x), -1)
    return(ddiff)}

fit_SAD<-function(totAb=1e7, rich=50, simpson=40, int_lwr=1e-4,int_uppr=1e2, dstr="lnorm"){
    #check feasibility
    if(simpson>rich){return("ERROR: Hill-Simpson diversity cannot be greater than richness")}

    #check dstr makes sense
    ifelse( !(dstr %in% c("lnorm", "gamma")), return("ERROR: dstr must be either `lnorm` or `gamma`"),

        #generate SAD when dstr=lnorm
        {if(dstr=="lnorm"){
            #uses an optimizer called uniroot to find x when ur_lognorm(x)==0
                fit_par=tryCatch(uniroot(function(x){ur_lognorm(simpson=simpson, rich=rich,x)}, lower=int_lwr, upper=int_uppr), error=function(e) message("ERROR: test int_lwr and int_uppr in ur_ function, output must have opposite signs"))

                #make sure to return rel abundances!
                abus=tryCatch(divers_lnorm(rich, fit_par$root), error=function(e) {message("did not fit param")
                    return(rep(100, length(rich)))}
                )
                abus<-abus/sum(abus)

        }

            #generate SAD when dstr="gamma"
        if(dstr=="gamma"){

            #uses an optimizer called uniroot to find x when ur_gamma(x)==0
            fit_par=tryCatch(uniroot(function(x){ur_gamma(simpson=simpson, rich=rich, x)}, lower=int_lwr, upper=int_uppr), error=function(e) message("ERROR: test int_lwr and int_uppr in ur_ function, output must have opposite signs"))

            #make sure to return rel abundances!
            abus=tryCatch(divers_gamma(rich,  fit_par$root), error=function(e) {
                message("did not fit param")
                return(rep(100, length(rich)))
                }
            )
            abus<-abus/sum(abus)
        }

            #return Hill-Shannon also
        shannon=dfun(abus, 0)
        # if(sum(abus==0)>0) print("WARNING: you simulated species with 0 abundance")

        return(list("distribution_info"=c("distribution"=dstr, "fitted parameter"=fit_par$root)
                    , "community_info"=c("richness"=rich, "Hill-Shannon"=shannon, "Hill-Simpson"=simpson)
                    , "rel_abundances"=abus

                   ))
        }
    )
}

#test function
# fit_SAD(dstr = "lnorm") #works
# fit_SAD(dstr = "nonsense")  #gives custom error (though not as error message)
# fit_SAD(dstr = "gamma") #works
# fit_SAD(dstr = "lnorm", rich = 50, simpson = 90) #gives custom error (though not as error message)
# fit_SAD(dstr = "lnorm", rich = 50, simpson = 2) #works
# fit_SAD(dstr = "gamma", rich = 50, simpson = 2) #works
# fit_SAD(dstr = "gamma", rich = 10, simpson = 6)



