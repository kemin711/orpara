########## Final Testing #######################
## not sensitive to b value 0.5 0.7 not different
source('build_gene.R')
source('loss.R')

testModern <- function(ng) {
	nrt <- 5000
	b <- 0.8
	conserve.par <- c(0.33, 0.82)
	### 1. use modern like ancestor
	alpha <- 2.4232631
	lambda <- 0.3130604
	fmodern.genome <- createFungalAncestor(ng, alpha, lambda, 1, conserve.par)
	#par(mfrow=c(3,1))
	#gs <- flatten.gene.feat(fmodern.genome);
	#drawhist(gs, 0)
	# loss parameter
	# for testing
	stat.res <- loss(fmodern.genome, b, nrt, 64)

	nrt <- 400
	fmodern2.genome <- createFungalAncestor(ng, alpha, lambda, 1, conserve.par)
	stat2.res <- loss(fmodern2.genome, b, nrt, 64)
# for looking at trend
	par(mfrow=c(2,2), mex=1.0, oma=c(2, 2, 1, 1))
	plot(stat.res$ril, stat.res$IPG, xlab='Relative Intron Location', ylab='IPG')
	tmpstat <- subset(stat.res, round>0)
	plot(log(tmpstat$numrt), tmpstat$IPG)
# second try
	plot(stat2.res$ril, stat2.res$IPG, xlab='Relative Intron Location', ylab='IPG')
	tmpstat2 <- subset(stat2.res, round>0)
    ipg.ril <- lm(IPG ~ ril, data=tmpstat2)
    abline(ipg.ril)
    predict.value <- predict(ipg.ril, newdata=data.frame(ril=0.5))
    error <- c(error, predict.value - stat2.res[1,5])
    print(paste('error: ', error))
	plot(log(tmpstat2$numrt), tmpstat2$IPG)
    
}
ng <- 10000
testModern(ng)


# return (estimate - real)
testModernError <- function(ng) {
	b <- 0.8
	conserve.par <- c(0.33, 0.82)
	alpha <- 2.4232631
	lambda <- 0.3130604
	nrt <- 380
	fmodern2.genome <- createFungalAncestor(ng, alpha, lambda, 1, conserve.par)
	stat2.res <- lossNoDraw(fmodern2.genome, b, nrt, 64)
	tmpstat2 <- subset(stat2.res, round>0)
    ipg.ril <- lm(IPG ~ ril, data=tmpstat2)
    predict.value <- predict(ipg.ril, newdata=data.frame(ril=0.5))
    est_real <- predict.value - stat2.res[1,5]
    print(paste('error: ', est_real))
	return(est_real)
}

ng <- 10000
error <- numeric()
for (i in 1:10) {
    error <- c(error, testModernError(ng))
}
error
mean(error)
sd(error)


#lm.x <- lm(tmpstat$IPG ~  log(tmpstat$numrt))
#abline(lm.x)

### 2. use more ancient fungal genome of 10 EPG
# 2. By RT method about 10 EPG
# since exon length will be to too long
# I will scale the alpha by 0.65
ng <- 10000
alpha <- 2.298
lambda <- 0.256
conserve.par <- c(0.37, 0.9)
b <- 0.8
nrt <- 8000
fancient.genome <- createFungalAncestor(ng, alpha, lambda, 0.5905, conserve.par)
# you only need to draw the ancestral genome once
# for other simulations comment it out
#drawhistGenome(fancient.genome)
# simulation result:
# RIL: 0.500, Exlen: 133.4 nt, IPG: 9.0
# CDS length 1200
# for testing
stat.res <- loss(fancient.genome, b, nrt, 64)
par(mfrow=c(2,1), mex=1.0, mar=c(4, 5, 2,2))
plot(stat.res$ril, stat.res$IPG, xlab='Relative Intron Location',
    ylab='IPG', ylim=c(1, 10))
# showing first data point
tmpstat <- subset(stat.res, round>0)
#tmpstat0 <- subset(stat.res, round==0)
#tmpstat0$numrt <- 1
#tmpstat <- rbind(tmpstat0, tmpstat)
#plot(log(tmpstat$numrt), tmpstat$IPG)
# not showing first data point
plot(log(tmpstat$numrt), tmpstat$IPG, xlab='ln(Number of RT)',
    ylab='IPG', ylim=c(1,10))

# simulate the beginning more closely
ng <- 10000
alpha <- 2.298
lambda <- 0.256
b <- 0.8
nrt <- 1000
fancient2.genome <- createFungalAncestor(ng, alpha, lambda, 0.5905, conserve.par)
# you only need to draw the ancestral genome once
# for other simulations comment it out
#drawhistGenome(fancient.genome)
# simulation result:
# RIL: 0.500, Exlen: 133.4 nt, IPG: 9.0
# CDS length 1200
# for testing
stat2.res <- loss(fancient.genome, b, nrt, 64)


# plotting two simulation summary in one figure.
par(mfrow=c(2,2), mex=1.0, mar=c(4, 5, 2,2))
plot(stat.res$ril, stat.res$IPG, xlab='Relative Intron Location',
    ylab='IPG', ylim=c(1, 10))
tmpstat <- subset(stat.res, round>0)
plot(log(tmpstat$numrt), tmpstat$IPG, xlab='ln(Number of RT)',
    ylab='IPG', ylim=c(1,10))
# second round
plot(stat2.res$ril, stat2.res$IPG, xlab='Relative Intron Location',
    ylab='IPG', ylim=c(2, 10))
tmpstat2 <- subset(stat2.res, round>0)
plot(log(tmpstat2$numrt), tmpstat2$IPG, xlab='ln(Number of RT)',
    ylab='IPG', ylim=c(2,10))


#lm.x <- lm(tmpstat$IPG ~  log(tmpstat$numrt))
#abline(lm.x)

### 2b. use more ancient fungal genome of 16 EPG
# 2. By RT method about 10 EPG
# since exon length will be to too long
# I will scale the alpha by 0.446
ng <- 10000
# alpha/lambda=15 => 16 exons per genes
alpha <- 3.225
lambda <- 0.215
b <- 0.7
nrt <- 8000
conserve.par <- c(0.3, 0.8)
fancient.genome <- createFungalAncestor(ng, alpha, lambda, 0.4459607, conserve.par)
gs <- flatten.gene.feat(fancient.genome);
par(mfrow=c(3,1))
drawhist(gs)
## use the average gene length to normalize the 
# scaling factor
# simulation result:
# RIL: 0.500, Exlen: 133.4 nt, IPG: 9.0
# CDS length 1200
# for testing
stat.res <- loss(fancient.genome, b, nrt, 64)
par(mfrow=c(2,1), mex=1.1)
plot(stat.res$ril, stat.res$IPG)
tmpstat <- subset(stat.res, round>0)
tmpstat0 <- subset(stat.res, round==0)
tmpstat0$numrt <- 1
tmpstat <- rbind(tmpstat0, tmpstat)
plot(log(tmpstat$numrt), tmpstat$IPG, xlim=c(0, 11))
#lm.x <- lm(tmpstat$IPG ~  log(tmpstat$numrt))
#abline(lm.x)


# 3. ancestral gene
# this is not very realistic
exlen.param <- c(1.1408, 0.01755)
# mean= 65.00285 nt exon length
numin.param <- c(4.003, 0.223)
# mean = 17.95 exons per gene
consv.param <- c(0.2, 0.8)
ng <- 2000
b <- 0.5
nrt <- 800
eancient.genome <- testEukaryoteAncestor(ng, numin.param, 
	exlen.param, consv.param)

stat.res <- loss(eancient.genome, b, nrt, 64)
par(mfrow=c(2,1), mex=1.1)
plot(stat.res$ril, stat.res$IPG)
tmpstat <- subset(stat.res, round>0)
tmpstat0 <- subset(stat.res, round==0)
tmpstat0$numrt <- 1
tmpstat <- rbind(tmpstat0, tmpstat)
plot(log(tmpstat$numrt), tmpstat$IPG, xlim=c(0, 11))



 
#-- blow are old
ng <- 10000
alpha <- 2.4232631
lambda <- 0.3130604
conserve.frac <- 0.5
conserve.per.gene <- 0.5

### combined parameters
numintron.par <- c(alpha,lambda)
exlen.param.file <- "gamma_param_exlen.tab"
conserve.param <- c(conserve.frac, conserve.per.gene)

gmodel <- build.gene(ng, numintron.par, conserve.param, exlen.param.file);
#####
#b <- 0.4
#nrt <- 5000
#gmnew <- gmodel
#par(mfrow=c(4,4), mex=0.5)
# at certain fraction of the genes have one conserved intron
# say 10%, different organisms have different 
# values for this one

b <- 0.9
nrt <- 5000
gmnew <- gmodel

ni2loc1 <- main(gmnew,b,5000, 100);
ni2loc1 <- c(0.5,7.74,ni2loc1);
nl <- matrix(ni2loc1, byrow=TRUE, ncol=2)
par(mfrow=c(1,1))
plot(nl[,1], nl[,2])

#gmnew <- gmodel
#ni2loc2 <- main(gmnew,b,50000);

ni2loc.more <- c(ni2loc1, ni2loc2)
nl <- matrix(ni2loc.more, byrow=TRUE, ncol=2)

par(mfrow=c(1,1))
plot(nl[,1], nl[,2])



