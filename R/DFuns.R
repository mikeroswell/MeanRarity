###### Here are the key mean rarity functions

###############################################################
##############################################################

#transformation and back-transformation functions
#' @export
pfun=function(x, pow){
  if (pow==0) return(log(x))
  r <- sign(pow)*(x)^pow
  return(r)
}

#' @export
ipfun=function(x, pow){
  if (pow==0) return(exp(x))
  r <- sign(pow)*(x)^(1/pow)
  return(r)
}

#' @export
dfun<-function(ab, l){
  ab<-ab[ab!=0]
  rp <- ab/sum(ab)
  if(l==0) {return(exp(sum(rp*log(1/rp))))}
  return(sign(l)*ipfun(sign(l)*sum(rp*pfun(1/rp, l)),l))
}
