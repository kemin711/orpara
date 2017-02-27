# exon length for genes with 15 introns
# parameter is (4.183, 0.02389)
# 

x <- seq(0, 700, 1)
alpha <- 4.183
lambda <- 0.02389
plot(x, dgamma(x))



#estimate the best parameter for 
# number of exons per gene by looking
# use number of introns per gene because this is in the range
# of gamma distribution.
# the parameter spror1 1.923, 0.313

