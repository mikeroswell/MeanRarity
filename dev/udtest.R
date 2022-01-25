
library(shellpipes)

loadEnvironments()

ll <- 1e-10
lh <- 1
steps <- 20

lvals <- exp(seq(log(ll), log(lh), length.out=steps))

for (l in lvals){
	print("l")
	print(l)
	print(down(up(1.1, l), l))
	print(down(up(2, l), l))
	print(down(up(10, l), l))
}

un <- numeric()
us <- numeric()
for (l in lvals){
	un <- up(10, l)
	un <- up(10, l, t=0)
}


