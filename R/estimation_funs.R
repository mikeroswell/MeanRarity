# Hodegepodge of code from Anne Chao group and Dushoff/Roswell for estimating Hill diversity and its uncertainty.
# Much of this is both pirated from Chao and not implemented in any of our work

# R scripts for computing diversity (Hill numbers) profile using individual-based abundance data or sampling-unit-based incidence data.
# In all functions, param x is a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (for incidence data).
# For incidence data, the first entry of x must be the number of sampling units.
# In all functions, param q is the diversity order; the suggested range for q is [0, 3].
# If you use the scripts for publishing papers, please cite Chao and Jost 2015 MEE paper (Appendix S8).

###########################################
#MR note: commenting out the incidence stuff b/c easier to focus on one thing at a time.

#-----------------------
# The empirical profile
#-----------------------

#This function probably not necessary since we have dfun, but this has incidence stuff which we could possibly care about.

#' Hill(x, q, datatype) is a function of obtaining the empirical Hill numbers of order q based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Hill <- function(x,q,datatype = c("abundance","incidence")){
  if(datatype=="incidence"){x = x[-1]}
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

#-----------------------------------------
# The bootstrap method for obtaining s.e.
#-----------------------------------------
#' Bt_prob_abu(x) is a function of estimating the species probabilities in the bootstrap assemblage based on abundance data.
#' @param x a vector of species sample frequencies.
#' @return a numeric vector.


###########################
# MR note: I had misrepresented the way that unobserved species were added to augmented sample. Slightly more nuanced than just adding unobserved as singletons.
Bt_prob_abu = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  #compute the coverage here
  C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
  #use coverage to define a weighting for observed frequencies... this is lambda_hat in Chao et al. 2013 and 2014 appendices explaining this bootstrapping procedure. I haven't quite wrapped my head around this yet.
  #coverage deficit=p(next individual is a new species)=proportion of true community absent from sample
  # the denominator is the expected probability of observing all x if x/n=p
  W = max((1-C)/sum(x/n*(1-x/n)^n), 0, na.rm=T)

  #use that weighting here to get p.new for observed species in x
  p.new = x/n*(1-W*(1-x/n)^n)
  #then get number of species observed 0 times using chao1
  f0 = ceiling(ifelse(f2>0, (n-1)/n*f1^2/(2*f2),max((n-1)/n*f1*(f1-1)/2, 0, na.rm=T))) #edited to deal with case where f1, f2=0
  #assume that all unobserved have equal p, given by total coverage deficit divided by number of unobserved.
  p0 = max((1-C)/f0, 0, na.rm=T) #consider issues if p0==0
  #p.new includes estimated p's for the observed plus estimated p's for unobserved.
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

## Can we do a better job of making a fake_community.md?
Bt_prob_abu_fiddle = function(x){
	x = x[x>0]
	n = sum(x)
	p <- x/n

	## What is the smallest number of unobserved species consistent
	## with a postulated coverage and target value of pair probability?
	umin <- function(p, C, pp){
		gap <- pp-C^2*sum(p^2)
		if (gap<=0) stop("C too large in umin (Bt_prob_abu_fiddle)")
		return((1-C)^2/gap)
	}
	## What is the expected richness of a sample from a community?
	## or pseudo-community?
	expRich <- function(n, alpha, X=NULL, u=1){
		return(
			sum(1-(1-alpha)^n)
			+ ifelse(is.null(X), 0, u*(1-(1-X/u)^n))
		)
	}
	## What is the expected richness of a hypothetical pseudo-community
	## with u=umin (not an integer) and even distribution therein?
	## Optionally subtract a postulated richness (for uniroot)
	uminRich <- function(C, n, p, pp, r=0){
		u <- umin(p, C, pp)
		return(expRich(n, C*p, 1-C, u)-r)
	}

	## Our indirect coverage estimate is a coverage that produces a
	## pseudo-community with expected richness equal to observed richness
	indCov <- function(x){
		eps <- 1e-3
		r <- length(x)
		n <- sum(x)
		pp <- sum(x*(x-1))/(n*n-1)
		Cmax <- pp/sum((x/n)^2)
		print(r)
		print(uminRich(eps*Cmax, n=n, p=x/n, pp=pp))
		print(uminRich((1-eps)*Cmax, n=n, p=x/n, pp=pp))
		return(uniroot(uminRich, lower=eps*Cmax, upper=(1-eps)*Cmax
			, n=n, p=x/n, pp=pp, r=r
		))
	}
	return(indCov(x))
}

## Sampled instead of deterministic abundance resampling
## Various things tried, but does not help with conservative CIs
Bt_prob_abu_samp = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = rpois(1, sum(x==1))
  f2 = rpois(1, sum(x==2))
  #Coverage
  C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
  #use coverage to define a weighting for observed frequencies... this is lambda_hat in Chao et al. 2013 and 2014 appendices explaining this bootstrapping procedure.
  #coverage deficit=p(next individual is a new species)=proportion of true community absent from sample
  # the denominator is the expected probability of observing all x if x/n=p
  W = (1-C)/sum(x/n*(1-x/n)^n)

  #use that weighting here to get p.new for observed species in x
  p.new = x/n*(1-W*(1-x/n)^n)
  #then get number of species observed 0 times using chao1
  f0 = rpois(1, ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  #assume that all unobserved have equal p, given by total coverage deficit divided by number of unobserved.
  p0 = (1-C)/f0
  #p.new includes estimated p's for the observed plus estimated p's for unobserved.
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob_inc(x) is a function of estimating the species incidence probabilities in the bootstrap assemblage based on incidence data.
#' @param x a vector of incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @return a numeric vector.

#' Bt_prob_inc = function(x){
#'     n = x[1]
#'     x = x[-1]
#'     U = sum(x)
#'     f1 = sum(x==1)
#'     f2 = sum(x==2)
#'     A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
#'     C=1-f1/U*(1-A)
#'     W=U/n*(1-C)/sum(x/n*(1-x/n)^n)
#'
#'     p.new=x/n*(1-W*(1-x/n)^n)
#'     f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
#'     p0=U/n*(1-C)/f0
#'     p.new=c(p.new,rep(p0,f0))
#'     return(p.new)
#' }
#'
#' #' Bt_prob(x,datatype) combines the two functions Bt_prob_abu and Bt_prob_inc for a specified datatype.
#' #' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' #' @param datatype a character of data type,"abundance" or "incidence".
#' #' @return a numeric vector.
#'
#' Bt_prob = function(x,datatype = c("abundance","incidence")){
#'     datatype = match.arg(datatype,c("abundance","incidence"))
#'     if(datatype == "abundance"){
#'         prob = Bt_prob_abu(x)
#'     }else{
#'         prob = Bt_prob_inc(x)
#'     }
#'     return(prob)
#' }

#' Bootstrap.CI(x,q,B,datatype,conf) is a function of calculating the bootsrapping standard error based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param B an integer to specify the number of replications in the bootstrap procedure, B = 1000 is suggested for constructing confidence intervals;
#'  To save running time, use a smaller value (e.g. B = 200)..
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 3 matrices including respectively the difference between the average and lower confidence bound of the B bootstrap estimates,
#'  the difference between the upper confidence bound and the average of the B bootstrap estimates, and the bootstrap standard error of the diversity estimate.
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

# x<-subsam(usersguide, 10)
# q<-0
# datatype<-"abundance"
# B<-1000

Bootstrap.CI = function(x,q,B = 1000,datatype = c("abundance","incidence"),conf = 0.95){
    datatype = match.arg(datatype,c("abundance","incidence"))
    p.new = Bt_prob_abu(x)
    # p.new = Bt_prob(x,datatype)
    n = ifelse(datatype=="abundance",sum(x),x[1])
    # set.seed(456)
    if(datatype=="abundance"){
        data.bt = rmultinom(B,n,p.new)
    }
    # }else{
    #     data.bt = rbinom(length(p.new)*B,n,p.new)
    #     data.bt = matrix(data.bt,ncol=B)
    #     data.bt = rbind(rep(n,B),data.bt)
    # }
    #
    mle = apply(data.bt,2,function(x)SpadeR:::Hill(x,q,datatype))
    #making this all for just abundance
    pro = apply(data.bt,2,function(x)SpadeR:::Chao_Hill_abu(x,q))

    mle.mean = mean(mle) #rowMeans(mle)
    pro.mean = mean(pro) #rowMeans(pro)

    #confidence intervals just based on quantiles of bootstraped distribution

    #confidence intervals for Hill diversity of sample
    # LCI.mle =  -apply(mle,1,function(x)quantile(x,probs = (1-conf)/2)) + mle.mean
    # UCI.mle = apply(mle,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - mle.mean
    #
    # #confidence intervals for Chao-estimated Hill diversity
    # LCI.pro =  -apply(pro,1,function(x)quantile(x,probs = (1-conf)/2)) + pro.mean
    # UCI.pro = apply(pro,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - pro.mean
    #
    LCI.mle = -quantile(mle,probs = (1-conf)/2) + mle.mean
    UCI.mle = quantile(mle,probs = 1-(1-conf)/2) - mle.mean

    #confidence intervals for Chao-estimated Hill diversity
    LCI.pro =  -quantile(pro,probs = (1-conf)/2) + pro.mean
    UCI.pro = quantile(pro,probs = 1-(1-conf)/2) - pro.mean

    LCI = rbind(LCI.mle,LCI.pro)
    UCI = rbind(UCI.mle,UCI.pro)

    # sd.mle = apply(mle,1,sd)
    # sd.pro = apply(pro,1,function(x)sd(x,na.rm = T))
    # se = rbind(sd.mle,sd.pro)
    #consider making this an easier data structure where things aren't lists of lists.
    #return(list(LCI=LCI,UCI=UCI,se=se))
    return(c(LCI.mle=LCI.mle, LCI.pro=LCI.pro, UCI.mle=UCI.mle, UCI.pro=UCI.pro))
}


#########################################
# ttest
# myt<-function(x, true){pt((mean(x)-true)/(sd(x)/sqrt(length(x))), df=length(x-1))}

################################################################
#MR function to
# x<-1:5
# B<-10
# l<-1
# truediv<-5

#'  Assess Chao and Jost 2015 Hill diversity CI
#'
#'  Given empirical abundances \code{x}, use the Chao and Jost 2015 MEE method
#'  to estimate the sampling distribution of the Hill diversity (given exponent
#'  \code{l}), and see where known, true diversity \code{truediv} falls in that
#'  distribution
#'
#'  This function could be extended to look at e.g. incidence-based estimations.
#'
#'  @param x Numeric vector of integer species abundances in a sample
#'  @param B Scalar, number of replicat bootstrap draws
#'  @param l Scalar, exponent for scaling rarity in computing Hill diversity
#'  @param truediv Scalar, known true Hill diversity of the pool from which
#'    sample is drawn, for comparison to estimated sampling distribution
#'  @param conf Scalar, target coverage probability of estimated CI
#'
#'  @return data.frame with estimated p-value for true diversity, the diversity
#'    values of the upper and lower estimated cofidence limits, the asymptotic
#'    Hill diversity estimate, and the empirical diversity
#'
#'  @export
checkchao<-function(x
                    , B
                    , l
                    , truediv
                    , conf=0.95){ #, truemu_n
  n<-sum(x)
  #columns of this matrix are replicate boostraps
  data.bt = stats::rmultinom(B,n,Bt_prob_abu(x))
  #get estimator for each bootstrapped sample
  pro = apply(data.bt,2,function(boot)SpadeR:::Chao_Hill_abu(boot,1-l))
  #mean correction
  pro<-pro-mean(pro)+Chao_Hill_abu(x, 1-l)

  #break ties
  less<-sum(pro<truediv)/length(pro)
  more<-(length(pro)-sum(pro>truediv))/length(pro)
  p<-runif(1, min(less, more), max(less, more))

  lower<-max(pro[which(min_rank(pro)<=max(floor(B*(1-conf)/2),1))])
  upper<-min(pro[which(min_rank(-pro)<=max(floor(B*(1-conf)/2),1))])

  return(data.frame(p=p
                , lower=lower
                , upper=upper
                , truediv=truediv
                , "chaoest"=SpadeR:::Chao_Hill_abu(x, 1-l)
                , "obsD"=dfun(x,l)
            )

  )}
