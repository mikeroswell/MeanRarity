
### Convenience functions to sample community data, estimate empirical mean diversities
# library(iNEXT)
# library(tidyverse)
# library(EntropyEstimation)
requireNamespace("furrr")
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

