# test speed of SpadeR::Diversity vs. its components
library(tictoc)
library(SpadeR)


tic()
SpadeR::Diversity(1:25)
toc() # 0.54 seconds

tic()
SpadeR:::Chat.Ind(1:25, sum(1:25))
toc() # 0.001 seconds
