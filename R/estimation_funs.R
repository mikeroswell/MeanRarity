#this is code that combines Chao et al and other stuff and is copied in more than one repo; consider what happens if they diverge.

# R scripts for computing diversity (Hill numbers) profile using individual-based abundance data or sampling-unit-based incidence data.
# In all functions, param x is a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (for incidence data).
# For incidence data, the first entry of x must be the number of sampling units. 
# In all functions, param q is the diversity order; the suggested range for q is [0, 3].
# If you use the scripts for publishing papers, please cite Chao and Jost 2015 MEE paper (Appendix S8). 

###########################################
#MR note: commenting out the incidence stuff b/c easier to focus on one thing at a time.

#-----------------------------------------------
# Diversity profile estimator (abundance data)
#-----------------------------------------------
#' Chao_Hill_abu(x, q) is a function of obtaining estimators of Hill numbers of order q based on abundance data.
#' @param x a vector of species sample frequencies. 
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @return a numerical vector of diversity. 

Chao_Hill_abu = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#' #-----------------------------------------------
#' # Diversity profile estimator (incidence data)
#' #-----------------------------------------------
#' #' Chao_Hill_inc(x, q) is a function of obtaining estimators of Hill numbers of order q based on incidence data.
#' #' @param x a vector of species incidence-based sample frequencies. The first entry of x must be the number of sampling units.
#' #' @param q a numeric or a vector of diversity order.
#' #' @return a numerical vector.
#' 
#' Chao_Hill_inc = function(x,q){
#'     n = x[1]
#'     x = x[-1];x = x[x>0]
#'     U = sum(x)
#'     f1 = sum(x==1)
#'     f2 = sum(x==2)
#'     p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
#'     r <- 0:(n-1)
#'     Sub <- function(q){
#'         if(q==0){
#'             sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
#'         }
#'         else if(q==1){
#'             A <- sum(x/U*(digamma(n)-digamma(x)))
#'             B <- ifelse(f1==0|p1==1,0,f1/U*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
#'             exp(A+B)*U/n
#'         }else if(abs(q-round(q))==0){
#'             A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
#'             ifelse(A==0,NA,((n/U)^q*A)^(1/(1-q)))
#'         }else {
#'             sort.data = sort(unique(x))
#'             tab = table(x)
#'             term = sapply(sort.data,function(z){
#'                 k=0:(n-z)
#'                 sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
#'             })
#'             A = sum(tab*term)
#'             B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
#'             ((n/U)^q*(A+B))^(1/(1-q))
#'         }
#'     }
#'     sapply(q, Sub)
#' }
#' 
#' #' Chao_Hill(x, q,datatype) combines Chao_Hill_abu and Chao_Hill_inc given a specified datatype (either abundance data or incidence data).
#' #' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' #' @param q a numeric or a vector of diversity order.
#' #' @param datatype a character of data type,"abundance" or "incidence".
#' #' @return a numerical vector.
#' 
#' Chao_Hill = function(x,q,datatype = c("abundance","incidence")){
#'     datatype = match.arg(datatype,c("abundance","incidence"))
#'     if(datatype == "abundance"){
#'         est = Chao_Hill_abu(x,q)
#'     }else{
#'         est = Chao_Hill_inc(x,q)
#'     }
#'     return(est)
#' }

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
    mle = apply(data.bt,2,function(x)Hill(x,q,datatype))
    #making this all for just abundance
    pro = apply(data.bt,2,function(x)Chao_Hill_abu(x,q))
    
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

#------------------------
# Main function ChaoHill
#------------------------
#' ChaoHill(dat, datatype, from, to, interval, B, conf) is the function of calculating the empirical and the proposed diversity profile, 
#' their bootsrap standard errors and confidance intervals.
#' @param dat a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param from a numeric number of diversity order q (the start order of profile).
#' @param to a numeric number of diversity order q (the end order of profile).
#' @param interval a numeric number to specify each increment of q from the start to end order.
#' @param B an integer to specify the number of bootstrap replications, B = 1000 is suggested for constructing confidence intervals; 
#'  To save running time, use a smaller value (e.g. B = 200).
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 4 matrices including respectively diversity estimates, bootstrap standard errors, lower confidence bounds, and upper confidence bounds. 
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

ChaoHill <- function(dat, datatype=c("abundance", "incidence"), from=0, to=3, interval=0.1, B=1000, conf=0.95){ 
  datatype = match.arg(datatype,c("abundance","incidence"))
  # for real data estimation
  
  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)
  q <- seq(from, to, by=interval)
  
  #-------------
  #Estimation
  #-------------
  MLE=Hill(dat,q,datatype)
  
  qD_pro=Chao_Hill(dat,q,datatype)
  
  CI_bound = Bootstrap.CI(dat,q,B,datatype,conf)
  se = CI_bound$se
  #-------------------
  #Confidence interval
  #-------------------
  tab.est=data.frame(rbind(MLE,qD_pro))
  
  LCI <- tab.est - CI_bound$LCI
  UCI <- tab.est + CI_bound$UCI
  
  colnames(tab.est) <- colnames(se) <- colnames(LCI) <- colnames(UCI) <- paste("q = ", q, sep="")    
  rownames(tab.est) <- rownames(se) <- rownames(LCI) <- rownames(UCI) <- c("Observed", "Chao_2013")
  return(list(EST = tab.est, 
              SD = se,
              LCI = LCI,
              UCI = UCI))
  
}


#----------------------------
# Plot of confidence interval
#----------------------------
#' conf.reg(x_axis,LCL,UCL,...) is a function to plot the confidence region.
#' 
#' \code{conf.reg} uses polygon to draw a confidence band plot
#' 
#' @param x_axis a vector of diversity orders.
#' @param LCL a vector of lower confidence bounds.
#' @param UCL a vector of upper confidence bounds.
#' @param ... further arguments to be passed to \code{polygon}
#' @return a polygon plot

conf.reg=function(x_axis,LCL,UCL,...) {
  x.sort <- order(x_axis)
  x <- x_axis[x.sort]
  LCL <- LCL[x.sort]
  UCL <- UCL[x.sort]
  polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
}

##############################################################################
### iNEXT coverage function
##############################################################################
Chat.Ind <- function(x, m=sum(x)){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  A <- ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
  Sub <- function(m){
    out <- 1-f1/n*A
  }
  sapply(m, Sub)		
}
###############################################################
##############################################################
# our diversity functions

#transformation and back-transformation functions
pfun=function(x, pow){
  if (pow==0) return(log(x))
  r <- sign(pow)*(x)^pow
  return(r)
}

ipfun=function(x, pow){
  if (pow==0) return(exp(x))
  r <- sign(pow)*(x)^(1/pow)
  return(r)
}
dfun<-function(ab, l){
  ab<-ab[ab!=0]
  rp <- ab/sum(ab)
  if(l==0) {return(exp(sum(rp*log(1/rp))))}
  return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/rp, l)),l))
}

### MR addition: simple communities to test doubling property
# library(iNEXT)
library(tidyverse)
# library(EntropyEstimation)
#function to take an abundance vector and subsample to size
subsam<-function(ab_vec, size=sum(ab_vec)){
  inds<-unlist(lapply(1:length(ab_vec), function(x){
    rep(x, ab_vec[x])
  }))
  sam<-sample(inds, size=size, replace=FALSE)
  ss<-unlist(lapply(1:length(ab_vec), function(y){
    length(which(sam==y))
  }))
  return(ss)
}
#take a subset of a bunch of community vectors, each subset of equal size 
subcom<-function(com, size){
  t(apply(com, 1, function(x){
    subsam(x, size)}
  ))
}

# Nielsen<-function(x){
#   w<-sum(x)
#   p<-x/sum(x)
#   d<-(w-1)^2/(sum(p^2)*(w+1)*(w-2)+3-2)
#   return(d)
# }

#estimates hill with order q for each comm at each size in sequence
raref<-function(from, to, by, comm,n=1, q){
  # ifelse(para==T, {
  require(furrr)
  nc<-detectCores()-1
  plan(strategy=multiprocess, workers=nc)
  p<-future_map_dfr(1:n, function(z){map_dfr(lapply(seq(from, to, by), function(b){
    o1<-apply(subcom(comm, b),1, function(x){
      # if(q==2){
      #     est<-Nielsen(x)
      # }
      # else{
      mrest<-fsd(ab=x,l=-q+1)
      est<-Chao_Hill_abu(x, q=q)#}
      emp<-dfun(ab=x, l=-q+1)
      
      coverage<-Chat.Ind(x)
      out<-rbind(divest=est, zhangest=mrest, divemp=emp, coverage=coverage, size=rep(b, length(est)), q=rep(1-q, length(est)))
      return(out)
    })
    return(data.frame(comm=row.names(comm), divest=o1[1,], divzhang=o1[2,], divemp=o1[3,], coverage=o1[4,], inds=o1[5,], ell=o1[6,]))
  }), rbind)}) #},
  
  # p<-map_dfr(1:n, function(z){map_dfr(lapply(seq(from, to, by), function(b){
  #     o1<-apply(subcom(comm, b),1, function(x){
  #         # if(q==2){
  #         #     est<-Nielsen(x)
  #         # }
  #         # else{
  #         est<-Chao_Hill_abu(x, q=q)#}
  #         emp<-dfun(ab=x, l=-q+1)
  #         coverage<-Chat.Ind(x)
  #         out<-rbind(divest=est, divemp=emp, coverage=coverage, size=rep(b, length(est)), q=rep(1-q, length(est)))
  #         return(out)
  #     })
  #     return(data.frame(comm=row.names(comm), divest=o1[1,], divemp=o1[2,], coverage=o1[3,], inds=o1[4,], ell=o1[5,]))
  # }), rbind)})
  return(p)
}

#fudge simpson to test performance

fsr<-function(x){
  (x-0.8)/(sum(x)-0.8)
}


fsd<-function(ab, l){
  ab<-ab[ab!=0]
  rp <- ab/sum(ab)
  fs<-fsr(ab)
  if(l==0) {return(exp(sum(rp*log(1/fs))))}
  return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/fs, l)),l))
}

# JD's function: takes vector of sample abundances and returns Simpson's estimator of Simpson's Diversity
sApp <- function(samp){
  n <- sum(samp)
  return(1/(
    sum((samp/n)*((samp-1)/(n-1)))
  ))
}
#########################################
# ttest
# myt<-function(x, true){pt((mean(x)-true)/(sd(x)/sqrt(length(x))), df=length(x-1))}

################################################################
#MR function to take abundance vector and "l" and return the quantile of B bootstrap iterations of Chao technique that true value falls on. Could be extended to include sample Hill as Chao seems to have
# x<-1:5
# B<-10
# l<-1
# truediv<-5

checkchao<-function(x, B, l, truediv, conf=0.95){ #, truemu_n
  n<-sum(x)
  #columns of this matrix are replicate boostraps
  data.bt = rmultinom(B,n,Bt_prob_abu(x))
  #get estimator for each bootstrapped sample
  pro = apply(data.bt,2,function(boot)Chao_Hill_abu(boot,1-l))
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
                , "chaoest"=Chao_Hill_abu(x, 1-l)
                , "obsD"=dfun(x,l)
            )
         
  )}

################
#truemu computes the emprical average sample diveristy under sampling without replacement
truemu<-function(comm, size, reps, l,...){
    sam<-replicate(reps, subsam(comm, size))
    return(
        mean(
            apply(sam,2, function(x){dfun(x, l)})
        )
    )
}



sample_infinite<-function(commSAD, size){
    namevec<-sample(1:length(commSAD), size=size, prob=commSAD, replace=T)
    mysam <- unlist(lapply(1:length(commSAD), function(y){
        length(which(namevec==y))}))#substample the whole community with # individuals=size)
    return(mysam)}

truemu_inf<-function(comm, size, reps, l,...){ #comm is abundance vector; size, reps, l all constants
    sam<-replicate(reps, sample_infinite(comm, size=size))
    
    return(
        mean(
            apply(sam,2, function(x){dfun(x, l)})
        )
    )
}

