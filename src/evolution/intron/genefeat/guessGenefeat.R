# exon length for genes with 15 introns
# parameter is (4.183, 0.02389)
# 

-- this one should be ok
x <- seq(0, 200, 0.1)
alpha <- 1.1408
lambda <- 0.01755
plot(x, dgamma(x, alpha, lambda), type='l')
avg <- alpha/lambda
md <- (alpha-1)/lambda
avg 
65.00285
md 
8.022792
# for 1200 CDS, there should be 18.46073 exons per gene (EPG)

#estimate the best parameter for 
# number of exons per gene by looking
# use number of introns per gene because this is in the range
# of gamma distribution.
# the parameter spror1 1.923, 0.313


x <- seq(0, 50, 0.1)
alpha <- 1.923
lambda <- 0.313


alpha <- 4.003
lambda <- 0.223

alpha <- 4.003
lambda <- 0.218
plot(x, dgamma(x, alpha, lambda), type='l')
avg <- alpha/lambda
md <- (alpha-1)/lambda
avg 
md

# this is the final result for ancient gene pool
exlen.param <- c(1.1408, 0.01755)
numin.param <- c(4.003, 0.223)

