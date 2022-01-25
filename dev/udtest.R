
library(shellpipes)

loadEnvironments()

ll <- 1e-10
lh <- 1
steps <- 20

lvals <- -exp(seq(log(ll), log(lh), length.out=steps))

for (l in c(lvals, -lvals)){
	cat("\n")
	print(down(up(1.1, l), l))
	print(down(up(2, l), l))
	print(down(up(10, l), l))
}

quit()

un <- numeric()
us <- numeric()
for (l in lvals){
	un <- up(10, l)
}


