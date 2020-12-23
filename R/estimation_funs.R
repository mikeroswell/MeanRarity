# Hodegepodge of code from Anne Chao group and Dushoff/Roswell for estimating Hill diversity and its uncertainty.

# Bt_prob_abu = function(x){
#   x = x[x>0]
#   n = sum(x)
#   f1 = sum(x==1)
#   f2 = sum(x==2)
#   #compute the coverage here
#   C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
#   #use coverage to define a weighting for observed frequencies... this is lambda_hat in Chao et al. 2013 and 2014 appendices explaining this bootstrapping procedure. I haven't quite wrapped my head around this yet.
#   #coverage deficit=p(next individual is a new species)=proportion of true community absent from sample
#   # the denominator is the expected probability of observing all x if x/n=p
#   W = max((1-C)/sum(x/n*(1-x/n)^n), 0, na.rm=T)
#
#   #use that weighting here to get p.new for observed species in x
#   p.new = x/n*(1-W*(1-x/n)^n)
#   #then get number of species observed 0 times using chao1
#   f0 = ceiling(ifelse(f2>0, (n-1)/n*f1^2/(2*f2),max((n-1)/n*f1*(f1-1)/2, 0, na.rm=T))) #edited to deal with case where f1, f2=0
#   #assume that all unobserved have equal p, given by total coverage deficit divided by number of unobserved.
#   p0 = max((1-C)/f0, 0, na.rm=T) #consider issues if p0==0
#   #p.new includes estimated p's for the observed plus estimated p's for unobserved.
#   p.new=c(p.new,rep(p0,f0))
#   return(p.new)
# }

# ## Can we do a better job of making a fake_community.md?
# Bt_prob_abu_fiddle = function(x){
# 	x = x[x>0]
# 	n = sum(x)
# 	p <- x/n
#
# 	## What is the smallest number of unobserved species consistent
# 	## with a postulated coverage and target value of pair probability?
# 	umin <- function(p, C, pp){
# 		gap <- pp-C^2*sum(p^2)
# 		if (gap<=0) stop("C too large in umin (Bt_prob_abu_fiddle)")
# 		return((1-C)^2/gap)
# 	}
# 	## What is the expected richness of a sample from a community?
# 	## or pseudo-community?
# 	expRich <- function(n, alpha, X=NULL, u=1){
# 		return(
# 			sum(1-(1-alpha)^n)
# 			+ ifelse(is.null(X), 0, u*(1-(1-X/u)^n))
# 		)
# 	}
# 	## What is the expected richness of a hypothetical pseudo-community
# 	## with u=umin (not an integer) and even distribution therein?
# 	## Optionally subtract a postulated richness (for uniroot)
# 	uminRich <- function(C, n, p, pp, r=0){
# 		u <- umin(p, C, pp)
# 		return(expRich(n, C*p, 1-C, u)-r)
# 	}

# 	## Our indirect coverage estimate is a coverage that produces a
# 	## pseudo-community with expected richness equal to observed richness
# 	indCov <- function(x){
# 		eps <- 1e-3
# 		r <- length(x)
# 		n <- sum(x)
# 		pp <- sum(x*(x-1))/(n*n-1)
# 		Cmax <- pp/sum((x/n)^2)
# 		print(r)
# 		print(uminRich(eps*Cmax, n=n, p=x/n, pp=pp))
# 		print(uminRich((1-eps)*Cmax, n=n, p=x/n, pp=pp))
# 		return(uniroot(uminRich, lower=eps*Cmax, upper=(1-eps)*Cmax
# 			, n=n, p=x/n, pp=pp, r=r
# 		))
# 	}
# 	return(indCov(x))
# }

# ## Sampled instead of deterministic abundance resampling
# ## Various things tried, but does not help with conservative CIs
# Bt_prob_abu_samp = function(x){
#   x = x[x>0]
#   n = sum(x)
#   f1 = rpois(1, sum(x==1))
#   f2 = rpois(1, sum(x==2))
#   #Coverage
#   C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
#   #use coverage to define a weighting for observed frequencies... this is lambda_hat in Chao et al. 2013 and 2014 appendices explaining this bootstrapping procedure.
#   #coverage deficit=p(next individual is a new species)=proportion of true community absent from sample
#   # the denominator is the expected probability of observing all x if x/n=p
#   W = (1-C)/sum(x/n*(1-x/n)^n)
#
#   #use that weighting here to get p.new for observed species in x
#   p.new = x/n*(1-W*(1-x/n)^n)
#   #then get number of species observed 0 times using chao1
#   f0 = rpois(1, ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
#   #assume that all unobserved have equal p, given by total coverage deficit divided by number of unobserved.
#   p0 = (1-C)/f0
#   #p.new includes estimated p's for the observed plus estimated p's for unobserved.
#   p.new=c(p.new,rep(p0,f0))
#   return(p.new)
# }



################################################################

# x<-1:5
# B<-10
# l<-1
# truediv<-5

#' Assess Chao and Jost 2015 Hill diversity CI
#'
#' Given empirical abundances \code{x}, use the Chao and Jost 2015 MEE method to
#' estimate the sampling distribution of the Hill diversity (given exponent
#' \code{l}), and see where known, true diversity \code{truediv} falls in that
#' distribution
#'
#' This function could be extended to look at e.g. incidence-based estimations.
#'
#' @param x Numeric vector of integer species abundances in a sample
#' @param B Scalar, number of replicate bootstrap draws
#' @param l Scalar, exponent for scaling rarity in computing Hill diversity
#' @param truediv Scalar, known true Hill diversity of the pool from which
#'    sample is drawn, for comparison to estimated sampling distribution
#' @param conf Scalar, target coverage probability of estimated CI
#'
#' @return data.frame with estimated p-value for true diversity, the diversity
#'    values of the upper and lower estimated confidence limits, the asymptotic
#'    Hill diversity estimate, and the empirical diversity
#'
#' @noRd
checkchao <- function(x
                    , B
                    , l
                    , truediv
                    , conf = 0.95){ #, truemu_n
  n <- sum(x)
  #columns of this matrix are replicate boostraps
  data.bt = stats::rmultinom(B, n, Bt_prob_abu(x))
  #get estimator for each bootstrapped sample
  pro = apply(data.bt
              , 2
              , function(boot) Chao_Hill_abu(boot, l))
  #mean correction
  pro <- pro - mean(pro) + Chao_Hill_abu(x, l)

  #break ties
  less <- sum(pro < truediv) / length(pro)
  more <- (length(pro) - sum(pro > truediv)) / length(pro)
  p <- stats::runif(1, min(less, more), max(less, more))

  lower <- max(pro[which(dplyr::min_rank(pro) <= max(floor(B * (1 - conf) / 2), 1))])
  upper <- min(pro[which(dplyr::min_rank(-pro) <= max(floor(B * (1 - conf) / 2) , 1))])

  return(data.frame(p = p
                , lower = lower
                , upper = upper
                , truediv = truediv
                , "chaoest" = Chao_Hill_abu(x, l)
                , "obsD" = dfun(x, l)
            )

  )}


#' Approximate CI for sample Hill diversity
#'
#' Computes Chao and Jost 2015's suggested CI for empirical Hill diversity,
#' sampling from an infinite pool (i.e. with replacement).
#'
#' @template l_template
#' @param SAD List, output from function `fit_SAD`
#' @param size Scalar, integer number of individuals to sample.
#' @param B Scalar, integer number of bootstrap replicates.
#' @param treumun = Scalar, true expected diversity for sample of size = `size`.
#' @param conf Scalar, nominal confidence level.
#' @param ... Arguments passed to other functions
#'
#' @return A data frame with upper and lower confidence limits and the p-value
#'   for the true diversity, computed both with and without the mean correction
#'   implemented in Chao and Jost 2015 MEE.
#'
#' @references
#' \insertRef{Chao2015}{MeanRarity}
#'
#'
#' @noRd
obscp_inf <- function(l = l
                      , size = size
                      , SAD = SAD
                      , B = 2000
                      , truemun = truemun
                      , conf = 0.95
                      , ...){
  sam <- sample_infinite(SAD$rel_abundances, size = size)
  data.bt = stats::rmultinom(B
                      , size
                      , Bt_prob_abu(sam)) #this genenerates "bootstrapped" samples
  obs <- rarity(sam, l) # observed diversity
  pro = apply(data.bt, 2, function(boot) rarity(boot, l)) #sample diversity for bootstraps
  pro_mc <- pro - mean(pro) + obs

  less <- sum(pro_mc < truemun) / length(pro_mc)
  more <- (length(pro_mc) - sum(pro_mc > truemun))/length(pro_mc)
  p <- stats::runif(1, min(less, more), max(less, more)) #this is a way to simulate the theoretical distribution of p-values in spite of ties.

  lower <- max(pro_mc[which(dplyr::min_rank(pro_mc) <= max(floor(B * (1 - conf) / 2), 1))])
  upper <- min(pro_mc[which(dplyr::min_rank(-pro_mc) <= max(floor(B * (1 - conf)/2), 1))])

  less_no_mc <- sum(pro < truemun) / length(pro)
  more_no_mc <- (length(pro) - sum(pro > truemun)) / length(pro)
  p_no_mc <- stats::runif(1, min(less, more), max(less, more))

  lower_no_mc <- max(pro[which(dplyr::min_rank(pro) <= max(floor(B * (1 - conf) / 2) ,1))])
  upper_no_mc <- min(pro[which(dplyr::min_rank(-pro) <= max(floor(B * (1 - conf)/2) ,1))])

  return(data.frame("p" = p
                    , "p_no_mc" = p_no_mc
                    , upper = upper
                    , lower = lower
                    , upper_no_mc = upper_no_mc
                    , lower_no_mc = lower_no_mc
                    , "truemu" = truemun
                    , "obsD" = obs
                    , "l" = l
                    , "size" = size ))

}


#' Approximate CI for sample Hill diversity
#'
#' Computes Chao and Jost 2015's suggested CI for empirical Hill diversity,
#' sampling from a finite pool without replacement.
#'
#' @template l_template
#' @template ab_template
#' @param size Scalar, integer number of individuals to sample.
#' @param B Scalar, integer number of bootstrap replicates.
#' @param treumun = Scalar, true expected diversity for sample of size = `size`.
#' @param conf Scalar, nominal confidence level.
#' @param ... Arguments passed to other functions
#'
#' @return A data frame with upper and lower confidence limits and the p-value
#'   for the true diversity, computed both with and without the mean correction
#'   implemented in Chao and Jost 2015 MEE.
#'
#' @seealso \url{https://doi.org/10.1111/2041-210X.12349}
#'
#' @noRd
obscp_obs <- function(l = l
                      , size = size
                      , ab = ab
                      , B = 2000
                      , truemun = truemun
                      , conf = 0.95
                      , ...){
  sam <- sample_finite(ab, size = size)
  data.bt = stats::rmultinom(B
                      , size
                      , Bt_prob_abu(sam)) #this genenerates "bootstrapped" samples
  obs <- rarity(sam, l) # observed diversity
  pro = apply(data.bt, 2, function(boot) rarity(boot, l)) #sample diversity for bootstraps
  pro_mc <- pro - mean(pro) + obs

  less <- sum(pro_mc < truemun) / length(pro_mc)
  more <- (length(pro_mc) - sum(pro_mc > truemun))/length(pro_mc)
  p <- stats::runif(1, min(less, more), max(less, more)) #this is a way to simulate the theoretical distribution of p-values in spite of ties.

  lower <- max(pro_mc[which(dplyr::min_rank(pro_mc) <= max(floor(B * (1 - conf) / 2), 1))])
  upper <- min(pro_mc[which(dplyr::min_rank(-pro_mc) <= max(floor(B * (1 - conf)/2), 1))])

  less_no_mc <- sum(pro < truemun) / length(pro)
  more_no_mc <- (length(pro) - sum(pro > truemun)) / length(pro)
  p_no_mc <- stats::runif(1, min(less, more), max(less, more))

  lower_no_mc <- max(pro[which(dplyr::min_rank(pro) <= max(floor(B * (1 - conf) / 2) ,1))])
  upper_no_mc <- min(pro[which(dplyr::min_rank(-pro) <= max(floor(B * (1 - conf)/2) ,1))])

  return(data.frame("p" = p
                    , "p_no_mc" = p_no_mc
                    , upper = upper
                    , lower = lower
                    , upper_no_mc = upper_no_mc
                    , lower_no_mc = lower_no_mc
                    , "truemu" = truemun
                    , "obsD" = obs
                    , "l" = l
                    , "size" = size ))

}

#' God's unbiased estimator
#'
#' The mean rarity perspective suggests a novel approach to estimating true
#' diversity from a sample, using the fact that the sample mean is an unbiased
#' estimate of the population mean. Although how to estimate a species’ rarity
#' from sample data remains unsolved, we test the promise of this idea using
#' what we call God’s estimator: All-knowing, God can plug in the true
#' population parameters for each species. Nevertheless, maybe for fun, God
#' chooses to compute the sample mean, weighting species’ true rarities by
#' their sampled frequency.
#'
#'
#' @param freqs Non-negative numeric vector of observed species frequencies
#' @param true_p Non-negative numeric vector of true, population-level species
#'   frequencies
#' @template l_template
#' @examples
#' GUE(freqs =  1:10/sum(1:10), true_p =  1:10/sum(1:10), l = 1) # true richness 10
#' GUE(freqs =  rep(1:5,2)/sum(1:10), true_p =  1:10/sum(1:10), l = 1) # true richness 10
#' @export

GUE <- function(freqs, true_p, l) {
  refreq = freqs/sum(freqs)
  retrue = true_p/sum(true_p)
  wts = refreq[freqs*true_p > 0]
  rabs = retrue[freqs*true_p > 0]
  ifelse((l == 0)
         , exp(sum(wts * log(1/rabs)))
         , (sum(wts * (1/rabs)^l))^(1/l))
}

