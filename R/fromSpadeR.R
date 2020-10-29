#' Estimate sample coverage
#'
#' Computes Chao and Jost 2012's estimate of sample coverage, assuming that
#' individuals are independently and randomly sampled. This function computes
#' rarefaction and extrapolation estimates of coverage given an integer sample
#' size \code{m}. We modified the source code to default to observed sample size
#' for convenience, and so that parameter names matched analogous use elsewhere
#' in the `MeanRarity` package. This function is copied directly from
#' `SpadeR:::Chat.Ind()` from the R package SpadeR 0.1.1 by Anne Chao, K. H. Ma,
#' T. C. Hsieh and Chun-Huo Chiu.
#'
#' @importFrom Rdpack reprompt
#' @references \insertRef{Chao2012a}{MeanRarity}
#'
#' @return Scalar between 0 and 1, estimated sample coverage.
#'
#' @template ab_template
#' @param m Scalar, sample size at which to estimate coverage.
#'
#' @export
#'
#' @examples
#'
#' # generate sample
#' abs <- sample_infinite(fit_SAD(rich = 50, simpson = 20)[[3]], 150)
#'
#' # estimate coverage of sample
#' Chat.Ind(abs)
#'
#' # estimate coverage under rarefaction
#' Chat.Ind(abs, m = 100)

Chat.Ind = function(ab, m = sum(ab)){ #modification was to default m to observed sample size, params match package idioms

  x <- ab[ab > 0] #edited to use ab instead of x
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (f1 > 0 & f2 > 0) {
    a = (n - 1) * f1/((n - 1) * f1 + 2 * f2)
  }
  if (f1 > 1 & f2 == 0) {
    a = (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
  }
  if (f1 == 1 & f2 == 0) {
    a = 0
  }
  if (f1 == 0) {
    a = 0
  }
  Sub <- function(m) {
    if (m < n)
      out <- 1 - sum(x/n * exp(lchoose(n - x, m) - lchoose(n -
                                                             1, m)))
    if (m == n)
      out <- 1 - f1/n * a
    if (m > n)
      out <- 1 - f1/n * a^(m - n + 1)
    out
  }
  sapply(m, Sub)
}

#' Generate "augmented" sample for bootstrapping
#'
#' Estimates true species frequencies (including unobserved species). The
#' estimated species frequencies are bootstrapped to generate CI for diversity
#' estimates following the method of Chao and Jost 2015 MEE. We modified the
#' source code from `SpadeR` so that parameter names matched analogous use
#' elsewhere in the `MeanRarity` package. This function is copied directly
#' from`SpadeR:::Chat.Ind()` from the R package SpadeR 0.1.1 by Anne Chao, K. H.
#' Ma, T. C. Hsieh and Chun-Huo Chiu.
#'
#' @references
#' \insertRef{Chao2015}{MeanRarity}
#'
#' @template ab_template
#' @return numeric vector of Estimated species frequencies
#'
#' @export
#' @examples
#'
#' # a sample
#' ab <- c(1, 1, 1, 2, 2, 2, 3, 4, 5, 7, 9, 12, 13, 21, 25, 50)
#'
#' # sample frequencies
#' freq <- (ab/sum(ab))
#' freq[order(freq, decreasing = TRUE)]
#'
#' #estimated true frequencies
#' etf <- Bt_prob_abu(ab)
#' etf[order(etf , decreasing = TRUE)]
#' # adjustment evident from 4th species on, 2 new rare species added
#'
#' sum(etf) #rel. abundances still sum to 1
#'
Bt_prob_abu = function(ab){
  x = ab[ab > 0]
  n = sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  C = 1 - f1/n * ifelse(f2 > 0, (n - 1) * f1/((n - 1) * f1 +
                                                2 * f2), ifelse(f1 > 0, (n - 1) * (f1 - 1)/((n - 1) *
                                                                                              (f1 - 1) + 2), 0))
  W = (1 - C)/sum(x/n * (1 - x/n)^n)
  p.new = x/n * (1 - W * (1 - x/n)^n)
  f0 = ceiling( ifelse(f2 > 0
                       , (n - 1)/n * f1^2/(2 * f2)
                       , max((n - 1)/n*f1*(f1 - 1)/2
                             , 0
                             , na.rm = T))) #edited to deal with case where f1, f2=0
  p0 = (1 - C)/f0
  p.new = c(p.new, rep(p0, f0))
  return(p.new)
}

#' Estimate asymptotic Hill diversity
#'
#' Computes an estimate of the true community diversity based on a finite sample
#' assuming that individuals are randomly and independently sampled, after Chao
#' and Jost 2015. `SpadeR` source code modified to match our use of \eqn{\ell}
#' instead of `q` to parameterize Hill diversity. Code taken from
#' `SpadeR:::Chao_Hill_abu()`From the R package SpadeR 0.1.1 by Anne Chao, K. H.
#' Ma, T. C. Hsieh and Chun-Huo Chiu
#'
#' @references
#' \insertRef{Chao2015}{MeanRarity}
#'
#' @return Scalar estimate of Hill diversity given \eqn{\ell}
#'
#' @template ab_template
#' @template l_template
#'
#' @export
#'
#' @examples
#'
#' # generate sample
#' abs <- sample_infinite(fit_SAD(rich = 50, simpson = 20)[[3]], 150)
#'
#' # estimate true Hill diversity
#' Chao_Hill_abu(abs, l = 1) # Chao1 estimate of richness lower bound
#' Chao_Hill_abu(abs, l = 0) # Asymptotic estimate of Hill-Shannon diversity
#' Chao_Hill_abu(abs, l = -1) # Asymptotic estimate of Hill-Simpson diversity
#'
Chao_Hill_abu = function(ab, l){ #modified param names according to package idiom
  q = 1 - l # modification from source
  x = ab[ab > 0] # modification from source
  n = sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  p1 = ifelse(f2 > 0, 2 * f2/((n - 1) * f1 + 2 * f2), ifelse(f1 >
                                                               0, 2/((n - 1) * (f1 - 1) + 2), 1))
  Sub <- function(q) {
    if (q == 0) {
      sum(x > 0) + (n - 1)/n * ifelse(f2 > 0, f1^2/2/f2,
                                      f1 * (f1 - 1)/2)
    }
    else if (q == 1) {
      r <- 1:(n - 1)
      A <- sum(x/n * (digamma(n) - digamma(x)))
      B <- ifelse(f1 == 0 | p1 == 1, 0, f1/n * (1 - p1)^(1 -
                                                           n) * (-log(p1) - sum((1 - p1)^r/r)))
      exp(A + B)
    }
    else if (abs(q - round(q)) == 0) {
      A <- sum(exp(lchoose(x, q) - lchoose(n, q)))
      ifelse(A == 0, NA, A^(1/(1 - q)))
    }
    else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data, function(z) {
        k = 0:(n - z)
        sum(choose(k - q, k) * exp(lchoose(n - k - 1,
                                           z - 1) - lchoose(n, z)))
      })
      r <- 0:(n - 1)
      A = sum(tab * term)
      B = ifelse(f1 == 0 | p1 == 1, 0, f1/n * (1 - p1)^(1 -
                                                          n) * (p1^(q - 1) - sum(choose(q - 1, r) * (p1 -
                                                                                                       1)^r)))
      (A + B)^(1/(1 - q))
    }
  }
  sapply(q, Sub)
}






#' Approximate CI for observed and asymptotic Hill diversity
#'
#' Functionally, a wrapper for `SpadeR:::Bootstrap_CI`, slimmed down to deal
#' with only abundance data and set to return a data.frame. This is the
#' approximate CI suggested by Chao and Jost 2015 MEE. Source code copied and
#' pasted from `SpadeR:::Bootstrap_CI` from the R package SpadeR 0.1.1 by Anne
#' Chao, K. H. Ma, T. C. Hsieh and Chun-Huo Chiu
#'
#' @references
#' \insertRef{Chao2015}{MeanRarity}


#' @param x a vector of species sample frequencies (for abundance data) or
#'   incidence-based sample frequencies (1st entry must be the number of
#'   sampling unit).
#' @template l_template
#' @param B an integer to specify the number of replications in the bootstrap
#'   procedure, B = 1000 is suggested for constructing confidence intervals;
#'  To save running time, use a smaller value (e.g. B = 200)..
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 3 matrices including respectively the
#'   difference between the average and lower confidence bound of the B
#'   bootstrap estimates, the difference between the upper confidence bound and
#'   the average of the B bootstrap estimates, and the bootstrap standard error
#'   of the diversity estimate. In each matrix, the first row gives the results
#'   for the empirical diversity, and the second row gives the results for the
#'   proposed diversity estimates. Columns give the results for different orders
#'   of q.
#'
#' @noRd

Bootstrap.CI_df = function(x
                           , l
                           , B = 1000
                           # , datatype = c("abundance","incidence")
                           , conf = 0.95){
  # datatype = match.arg(datatype, c("abundance","incidence"))
  p.new = Bt_prob_abu(x)
  # p.new = Bt_prob(x,datatype)
  n = sum(x) # removed ifelse(datatype == "abundance",
  # set.seed(456)
  # if(datatype == "abundance"){
    data.bt = stats::rmultinom(B, n, p.new)
  # }
  # }else{
  #     data.bt = rbinom(length(p.new)*B,n,p.new)
  #     data.bt = matrix(data.bt,ncol=B)
  #     data.bt = rbind(rep(n,B),data.bt)
  # }
  #
  mle = apply(data.bt, 2, function(x)rarity(x, l)) #our Hill diversity function
  #making this all for just abundance
  pro = apply(data.bt, 2, function(x)Chao_Hill_abu(x, l)) #modified SpadeR function

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
  LCI.mle = -stats::quantile(mle, probs = (1 - conf) / 2) + mle.mean
  UCI.mle = stats::quantile(mle, probs = 1 - (1 - conf) / 2) - mle.mean

  #confidence intervals for Chao-estimated Hill diversity
  LCI.pro =  -stats::quantile(pro, probs = (1 - conf) / 2) + pro.mean
  UCI.pro = stats::quantile(pro, probs = 1 - (1 - conf) / 2) - pro.mean

  LCI = rbind(LCI.mle, LCI.pro)
  UCI = rbind(UCI.mle, UCI.pro)

  # sd.mle = apply(mle,1,sd)
  # sd.pro = apply(pro,1,function(x)sd(x,na.rm = T))
  # se = rbind(sd.mle,sd.pro)
  #consider making this an easier data structure where things aren't lists of lists.
  #return(list(LCI=LCI,UCI=UCI,se=se))
  return(c(LCI.mle = LCI.mle, LCI.pro = LCI.pro, UCI.mle = UCI.mle, UCI.pro = UCI.pro))
}

