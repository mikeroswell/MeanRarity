
library(shellpipes)

## Numerically stable inverse functions for calculating scaled exponents.
## The idea is that we up before down, so the concern is that ell might be small

## Transform so that u ~ (x^ℓ-1)/ℓ

up <- function(x, l, t=1e-4, t1=NULL){
	if(is.null(t1)) t1=t^2
	if(abs(l) <= t1) return(log(x))
	if(abs(l) < t) return(log(x) + l*log(x)^2/2)
	return((x^l-1)/l)
}

down <- function(u, l, t=1e-4, t1=NULL){
	if(is.null(t1)) t1=t^2
	if(abs(l) <= t1) return(exp(u))
	if(abs(l) < t) return(exp((sqrt(1+2*l*u)-1)/l))
	return((l*u+1)^(1/l))
}


pfun = function(x, pow){
  if (pow == 0) return(log(x))
  r <- sign(pow) * up(x,pow)
  return(r)
}

ipfun = function(x, pow){
  if (pow == 0) return(exp(x))
  x <- ifelse(sign(pow) * x < 0, 0, x) #added so that ggplot padding doesn't introduce negative values to scale
  r <- (sign(pow) * down(x, pow))
  return(r)
}

# modified rarity for better numerical behavior (precise rarity)
prec_rarity = function(ab, l, q = NULL){
  if(!is.null(q)){
    l = 1-q
    warning("l has been set to 1-q")
  }
  ab = ab[ab != 0]
  rp = ab/sum(ab)
  if(l == 0){return(exp(sum(rp * log(1/rp))))}
  return(down(sum(rp * up(1/rp, l)), l)) # replaces pfun and ipfun with up, down
}

# test that it gives finite answers
prec_rarity(c(10,1), l = 296)
prec_rarity(c(10,1), l = 297)
