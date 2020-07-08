
### Convenience functions to sample community data, estimate empirical mean diversities


#' Take an abundance vector and subsample to size
#'
#' Take a finite sample of individuals without replacement from a finite abundance vector.
#'
#' @param ab_vec A numeric vector of species abundances
#' @param size Number of individuals to sample, defaults to all of them, a scalar
#'
#' @return A numeric vector of species abundances, including 0's
#'
#' @export
#' @examples subsam(1:9, 15)
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


#' Subsample of several community vectors

#' A wrapper of \code{\link{subsam}} to take a subset of a bunch of communities,
#' each subset of equal abundance.
#'
#' @param com A list of abundance vectors.
#' @param size Number of individuals to sample, a scalar.
#'
#' @export
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


#' Estimate Hill diversity with order l=1-q under rarefaction
#'
#' This is a function (currently run in parallel with
#'      \code{\link[parallel]{detectCores}} and \code{\link[furrr]{future_map_dfr}} )
#'      that returns rarefied Hill diversity estimates for a list of sample (or true)
#'      abundance vectors.
#'
#'
#' Note to developer:  this might be the only place parallel is used, figure out
#'      if necessary and then either allow to not always be parallelized or consider
#'      omiting entirely
#'
#' @param from Scalar, smallest sample size in rarefaction
#' @param to Scalar, largest sample size in rarefaction
#' @param by Scalar, increment in \code{seq(from, to, by)}
#' @param comm
#' @param n
#' @param l exponent for scaling mean rarity. Scalar.
#' @param cores optional argument to set number of cores for parallel computing,
#'   defaults to \code{parallel::detectCores()-1}
#'
#' @return data.frame with various Hill-Diversity estimates and sample coverage
#'   estimates for each sample size in rarefaction
#'
#' @seealso \code{\link{Chat.Ind}}, \code{\link{subcom}}, \code{\link{dfun}},
#'   \code{\link{Chao_Hill_abu}},
#'
#' @noRd
#'

raref<-function(from, to, by, comm, n = 1, l, cores = NULL){
  # ifelse(para==T, {
  nc<-parallel::detectCores()-1
  furrr::plan(strategy=multiprocess, workers=ifelse(is.null(cores), nc, cores))
  p<-furrr::future_map_dfr(1:n, function(z){
    purrr::map_dfr(lapply(seq(from, to, by), function(b){
      o1<-apply(subcom(comm, b),1, function(x){
      # if(q==2){
      #     est<-Nielsen(x)
      # }
      # else{
      mrest<-fsd(ab=x,l=l)
      est<-SpadeR::Chao_Hill_abu(x, q=1-l)#}
      emp<-dfun(ab=x, l=l)

      coverage<-iNEXT:::Chat.Ind(x)
      out<-rbind(divest=est, zhangest=mrest, divemp=emp, coverage=coverage, size=rep(b, length(est)), q=rep(1-q, length(est)))
      return(out)
    })
    return(data.frame(comm=row.names(comm), divest=o1[1,], divzhang=o1[2,], divemp=o1[3,], coverage=o1[4,], inds=o1[5,], ell=o1[6,]))
  }), rbind)}) #},

  # p<-purrr:map_dfr(1:n, function(z){map_dfr(lapply(seq(from, to, by), function(b){
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



################
#' Compute emprical average sample diversity
#'
#' Based on replicate samples from a finite pool (samples taken without replacment).
#'
#' @param ab Numeric vector of species abundances.
#' @param size Scalar, number of individuals in sample
#' @param reps Scalar, number of replicate samples to take
#' @param l Scalar, exponent determining type of mean rarity
#'
#' @return scalar, emprical measure of the mean sample diversity from a larger pool
#'
#' @seealso \code{\link{dfun}}; \code{\link{subsam}}
#'
#' @noRd
truemu<-function(ab, size, reps, l,...){
  sam<-replicate(reps, subsam(ab, size))
  return(
    mean(
      apply(sam,2, function(x){dfun(x, l)})
    )
  )
}


#'  Sample a [relative] abundance vector with replacement
#'
#'  @param commSAD Numeric vector of [relative] abundances
#'  @param size Scalar, number of individuals in finite samples
#'
#'  @return A vector of integer species abundances
#'
#'  @export
sample_infinite<-function(commSAD, size){
  namevec<-sample(1:length(commSAD), size=size, prob=commSAD, replace=T)
  mysam <- unlist(lapply(1:length(commSAD), function(y){
    length(which(namevec==y))}))#substample the whole community with # individuals=size)
  return(mysam)}

#' Estimate empirical mean rarity (given l) for a finite-sized sample
#'
#' Samples taken from a species abundance distribution (sampling with replacement)
#'
#' @param comm numeric vector of species [relative] abundances
#' @param size scalar, number of individuals in each sample
#' @param reps Scalar, number of replicate samples
#' @param l Scalar, exponent for scaling mean rarity
#'
#' @return Scalar, empirical mean sample diversity given sampling with replacement
#'
#' @noRd
truemu_inf<-function(comm, size, reps, l,...){ #comm is abundance vector; size, reps, l all constants
  sam<-replicate(reps, sample_infinite(comm, size=size))
  return(
    mean(
      apply(sam,2, function(x){dfun(x, l)})
    )
  )
}

