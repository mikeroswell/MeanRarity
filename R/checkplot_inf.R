## checkplot_inf.R - compiled by RoxygenReady, a package by @vertesy


#' Generate the data for checkplots for Hill diversity
#'
#'  Estimating sampling uncertainty for Hill diversity estimates is an unsolved
#'  problem, but \insertCite{Chao2015}{MeanRarity} provide a heuristic solution. This function
#'  can be applied to a set species abundance distribution (\code{SAD}) to
#'  produce true and estimated Hill diversity sampling uncertainty (assuming
#'  individuals are randomly and independently sampled).
#'
#'
#' @param SAD List, output from MeanRarity::fit_SAD
#' @param B Integer number of replicate samples to take when bootstrapping for
#' the estimated uncertainty
#' @inherit ell_template
#' @param inds Integer number of individuals to sample from the SAD (with
#' replacement)
#' @param reps Integer number of replicate samples to take to generate sampling
#' uncertainty and checkplots/slugplots
#'
#' @references
#' \insertAllCited{}
#'
#' @noRd
checkplot_inf <-function (SAD, B = 2000, l, inds, reps) {
	hillname <- ifelse(l == -1, "Hill-Simpson", ifelse(l == 0, "Hill-Shannon", "richness"))
	td <- SAD$community_info[hillname]
	furrr::future_map_dfr(1:reps, function(x) {
		obs <- sample_infinite(SAD$rel_abundances, size = inds)
		chaotile <- checkchao(x = obs, B = B, l = l, truediv = td)
		return(myout = data.frame(p = chaotile$p, truediv = chaotile$truediv, chaoest = chaotile$chaoest,
			obsD = chaotile$obsD, upper = chaotile$upper, lower = chaotile$lower, l, inds, reps))
	})
}


