
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

saveEnvironment()
