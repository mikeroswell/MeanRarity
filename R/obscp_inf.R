## obscp_inf.R - compiled by RoxygenReady, a package by @vertesy


#' Observed diversity confidence interval
#'
#' Compute confidence interval for the expected diversity of a fixed sample,
#' taken with replacement
#'
#'
#'
#'
#'
#'
#' @inherit ell_template
#' @template q_template
#' @param size Integer number of individuals to sample from the SAD (with
#' replacement)
#' #' @param SAD List, output from MeanRarity::fit_SAD
#' @param B Integer number of replicate samples to take when bootstrapping for
#' the estimated uncertainty
#' @param truemun Scalar, known expected sample diversity
#' @param conf Scalar between 0 and 1, nominal confidence level.
#' @param ...
#'
#' @return Dataframe containing various estimated p-values, confidence limits,
#' sample and community characteristics
#'
#' @noRd

obscp_inf <- function(l = l
                      , q = NULL
                      , size = size
                      , SAD = SAD
                      , B = 2000
                      , truemun = truemun
                      , conf = 0.95, ...){
  if(!is.null(q)){
    l = 1-q
    warning("l has been set to 1-q")
  }
	sam = sample_infinite(SAD$rel_abundances, size = size)
	data.bt = stats::rmultinom(B, size, Bt_prob_abu(sam))
	obs = dfun(sam, l)
	pro = apply(data.bt, 2, function(boot) dfun(boot, l))
	pro_mc = pro - mean(pro) + obs
	less <- sum(pro_mc < truemun)/length(pro_mc)
	more = (length(pro_mc) - sum(pro_mc > truemun))/length(pro_mc)
	p = stats::runif(1, min(less, more), max(less, more))
	lower = max(pro_mc[which(dplyr::min_rank(pro_mc) <= max(floor(B * (1 - conf)/2), 1))])
	upper = min(pro_mc[which(dplyr::min_rank(-pro_mc) <= max(floor(B * (1 - conf)/2), 1))])
	less_no_mc = sum(pro < truemun)/length(pro)
	more_no_mc = (length(pro) - sum(pro > truemun))/length(pro)
	p_no_mc = stats::runif(1, min(less, more), max(less, more))
	lower_no_mc = max(pro[which(dplyr::min_rank(pro) <= max(floor(B * (1 - conf)/2), 1))])
	upper_no_mc = min(pro[which(dplyr::min_rank(-pro) <= max(floor(B * (1 - conf)/2), 1))])
	return(data.frame(p = p, p_no_mc = p_no_mc, upper = upper, lower = lower, upper_no_mc = upper_no_mc,
		lower_no_mc = lower_no_mc, truemu = truemun, obsD = obs, l = l, size = size))
}

