#' Take an abundance vector and subsample to size
#'
#' Take a finite sample of individuals without replacement from an integer
#' abundance vector.
#'
#' @template ab_template
#' @param size Number of individuals to sample, defaults to all of them, a
#'   scalar
#'
#' @return A numeric vector of species abundances, including 0's
#'
#' @seealso \code{\link{sample_infinite}}
#'
#' @concept Sampling
#'
#' @export
#' @examples sample_finite(1:9, 15)
sample_finite <- function(ab, size = sum(ab)){
  inds <- unlist(lapply(1:length(ab), function(x){
    rep(x, ab[x])
  }))
  sam <- sample(inds, size = size, replace = FALSE)
  ss <- unlist(lapply(1:length(ab), function(y){
    length(which(sam == y))
  }))
  return(ss)
}

subsam <- sample_finite


#' Subsample of several community vectors

#' A wrapper of \code{\link{sample_finite}} to take a subset of a bunch of
#' communities, each subset of equal abundance.
#'
#' @template ab_template
#' @param size Number of individuals to sample, a scalar.
#'
#' @concept Sampling
#'
#' @export
subcom <- function(ab, size){
  t(apply(ab, 1, function(x){
    sample_finite(x, size)}
  ))
}

# Nielsen<-function(x){
#   w<-sum(x)
#   p<-x/sum(x)
#   d<-(w-1)^2/(sum(p^2)*(w+1)*(w-2)+3-2)
#   return(d)
# }


#' Estimate Hill diversity with order \code{l = 1-q} under rarefaction
#'
#' This is a function (currently run in parallel with
#' \code{\link[parallel]{detectCores}} and \code{\link[furrr]{future_map_dfr}})
#' that returns rarefied Hill diversity estimates for a list of sample (or true)
#' abundance vectors.
#'
#'
#' Note to developer:  this might be the only place parallel is used, figure out
#' if necessary and then either allow to not always be parallelized or consider
#' omitting entirely
#'
#' @param from Scalar, smallest sample size in rarefaction.
#' @param to Scalar, largest sample size in rarefaction.
#' @param by Scalar, increment in \code{seq(from, to, by)}.
#' @param comm List of integer abundance vectors.
#' @param n Integer number of replicate rarefaction samples.
#' @template l_template
#' @param cores optional argument to set number of cores for parallel computing,
#'   defaults to \code{parallel::detectCores()-1}.
#' @param ... Additional arguments passed to other functions.
#'
#' @return data.frame with various Hill-Diversity estimates and sample coverage
#'   estimates for each sample size in rarefaction
#'
#' @seealso \code{\link{Chat.Ind}}, \code{\link{subcom}}, \code{\link{dfun}},
#'   \code{\link{Chao_Hill_abu}},
#'
#' @noRd
#'

raref <- function(from, to, by, comm, n = 1, l, cores = NULL){
  # ifelse(para==T, {
  nc <- parallel::detectCores() - 1
  future::plan(strategy = future::multiprocess
               , workers = ifelse(is.null(cores), nc, cores))
  p <- furrr::future_map_dfr(1:n, function(z){
    purrr::map_dfr(lapply(seq(from, to, by), function(b){
      o1 <- apply(subcom(comm, b), 1, function(x){
        mrest <- fsd(ab = x, l = l)
        est <- Chao_Hill_abu(x, l = l)#}
        emp <- rarity(ab = x, l = l)
        coverage <- Chat.Ind(x)
        out <- rbind(divest = est
                     , zhangest = mrest
                     , divemp = emp
                     , coverage = coverage
                     , size = rep(b, length(est))
                     , l = rep(l, length(est)))
        return(out)
    })
    return(data.frame(
      comm = row.names(comm)
      , divest = o1[1,]
      , divzhang = o1[2,]
      , divemp = o1[3,]
      , coverage = o1[4,]
      , inds = o1[5,]
      , ell = o1[6,]))
  }), rbind)})
  return(p)
}



################
#' Compute empirical average sample diversity
#'
#' Based on replicate samples from a finite pool (samples taken without
#' replacement).
#'
#' \code{ab} must be an integer vector.
#'
#' @template ab_template
#' @param size Scalar, number of individuals in sample
#' @param reps Scalar, number of replicate samples to take
#' @template l_template
#' @param ... Additional arguments passed to other functions.
#'
#' @return scalar, empirical measure of the mean sample diversity from a larger
#'   pool
#'
#' @seealso \code{\link{MeanRarity}}; \code{\link{subsam}}
#'
#' @noRd
truemu <- function(ab, size, reps, l, ...){
  sam <- replicate(reps, subsam(ab, size))
  return(
    mean(
      apply(sam, 2, function(x){dfun(x, l)})
    )
  )
}


#' Sample a \[relative\] abundance vector with replacement
#'
#' Subsample the whole community with number of individuals = \code{size}.
#'
#' @template ab_template
#' @param size Scalar, number of individuals in finite samples.
#'
#' @return A vector of integer species abundances
#'
#' @seealso \code{\link{sample_finite}}
#'
#' @concept Sampling
#'
#' @export
sample_infinite <- function(ab, size){
  namevec <- sample(1:length(ab)
                    , size = size
                    , prob = ab
                    , replace = TRUE)
  mysam <- unlist(lapply(1:length(ab), function(y){
    length(which(namevec == y))}))
  return(mysam)}

#' Estimate empirical mean rarity (given l) for a finite-sized sample.
#'
#' Samples taken from a species abundance distribution (sampling with
#' replacement).
#'
#' \code{ab} Can be absolute or relative abundances (need not be integers).
#'
#' @template ab_template
#' @param size scalar, number of individuals in each sample
#' @param reps Scalar, number of replicate samples
#' @template l_template
#'
#' @return Scalar, empirical mean sample diversity given sampling with
#'   replacement
#'
#' @noRd
truemu_inf <- function(ab, size, reps, l, ...){ #ab is abundance vector; size, reps, l all constants
  sam <- replicate(reps, sample_infinite(ab, size = size))
  return(
    mean(
      apply(sam,2, function(x){dfun(x, l)})
    )
  )
}

