-- for Sporo1 whose parameter is simulated by gamma
-- distribution: alpha=1.9232631, lambda <- 0.3130604
-- This is the number of exons distribution
-- we assume this as the ancestral state
-- assume 10,000 genes in the genome

# number of introns per gene
# the parameters for Sporo1, which has on average, 6.2 introns
# The conserved gens of Lacbi1 has 6.9 introns, so we we assume
# the ancestor has 7.8 introns
alpha <- 1.9232631
lambda <- 0.3130604
x <- c(seq(0,1, 0.02), seq(1.02, 70, 0.5))
plot(x, dgamma(x, alpha, lambda))
lines(x, dgamma(x, alpha+0.5, lambda))

### given a list of genes it 
# build the location into one single vector
# get the intron location and number of intron statistics
statgene <- function(gg) {
	N <- length(gg)
	loc <- numeric(); numin <- numeric();
	for (i in 1:N) {
		if (gg[[i]][1] > 0) {
			loc <- c(loc, gg[[i]])
			numin <- c(numin, length(gg[[i]]))
		}
		else {
			numin <- c(numin, 0)
		}
	}
	list(rel.loc=loc, num.intron=numin);
}

# returns the newly modified gene list
# will lose any selected intron 
intron.loss <- function(gg, b, nrt) {
	# pick a random list of genes to work on
	numgene <- length(gg);
	gi <- floor(runif(nrt, 1, numgene)+0.5)
	U <- runif(nrt)
	# loss random variable, from 5' to 3' loss pdf
	X <- 0.5*(sqrt(b^2+4*(1-b)*U)-b)/(1-b)
	i <- 1
	while (i <= nrt) {
		pos <- gi[i]
		stop.at <- X[i]
		if (gg[[pos]][1] != 0) { 
			tmp <- gg[[pos]];
			if (length(tmp[tmp < stop.at]) < 1) {
				gg[[pos]] <- c(0);
			}
			else {
				gg[[pos]] <- tmp[tmp<stop.at]
			}
		}
		i <- i+1
	}
	gg
}

### new  version that will use gmodel as input
# it preserved conserved introns
intron.loss2 <- function(gm, b, nrt) {
	numgene <- length(gm$reloc);
	gi <- floor(runif(nrt, 1, numgene)+0.5)
	U <- runif(nrt)
	X <- 0.5*(sqrt(b^2+4*(1-b)*U)-b)/(1-b)
	i <- 1
	while (i <= nrt) {
		pos <- gi[i]; stop.at <- X[i];
		print(paste('position:', pos))
		print(gm$reloc[[pos]])
		print(gm$reloc[[pos]][1])
		if (gm$reloc[[pos]][1] != 0) { 
			tmp <- gm$reloc[[pos]];
			savedidx <- tmp < stop.at;
			if (gm$conserve[i] > 0) { 
				print('Conserved savedidx')
				print(savedidx)
				savedidx[gm$conserve[i]]=T; 
			}
			print(savedidx)
			res <- tmp[savedidx];
			if (length(res) == 0) { gm$reloc[[pos]] <- list(0); }
			else { gm$reloc[[pos]] <- res; }
		}
		i <- i+1
	}
	gm
}

#######################################################################
## simulation of intron loss
ng <- 10000
consvfrac <- 0.2
#-- we assume the Sporolomyces lost 0.5 intron on average
# the number of introns per gene
numintron <- floor(rgamma(ng, alpha+0.5, lambda)+0.5)
genes <- list()
consvi <- numeric()
for ( i in 1:ng) {
	ci <- 0; v <- 0
	if (numintron[i] > 0) {
		v <- c(runif(numintron[i]))
		#l <- list(c(runif(numintron[i])))
		if (runif(1) < consvfrac) {
			ci <- floor(runif(1, 1, length(v))+0.5)
		}
	}
	genes <- c(genes, list(v)); consvi <- c(consvi, ci)
}
gmodel <- list(reloc=genes, conserve=consvi)


##### no conserved ntrons #######################
h <- histloc(genes)
hist(h, 100)

b <- 0.1
nrt <- 5000
genes.new <- genes
par(mfrow=c(4,4), mex=0.5)
# at certain fraction of the genes have one conserved intron
# say 10%, different organisms have different 
# values for this one

#for (i in 1:96) {
#for (i in 1:8) {
#for (i in 1:32) {
#for (i in 1:64) {
for (i in 1:199) {
	genes.new <- intron.loss(genes.new, b, nrt);
	if (i %% 8 == 0 & i*nrt > 320000) {
		h <- statgene(genes.new)
		hist(h$rel.loc, 100, main='', xlab='', ylab='', xlim=c(0,1), freq=F)
		hh <- hist(h$num.intron, 100, main='', xlab='', ylab='')
		avg.intron <- mean(h$num.intron)
		text(max(hh$breaks)*0.5, max(hh$count)*0.9, 
			labels=paste(i*nrt, "\n", sprintf("%1.2f",avg.intron)), 
			cex=0.9, adj=0)
	}
}

##### has conserved ntrons #######################
h <- histloc(genes)
hist(h, 100)

b <- 0.1
nrt <- 5000
gmnew <- gmodel
par(mfrow=c(4,4), mex=0.5)
# at certain fraction of the genes have one conserved intron
# say 10%, different organisms have different 
# values for this one

#for (i in 1:96) {
for (i in 1:8) {
#for (i in 1:32) {
#for (i in 1:64) {
#for (i in 1:199) {
	gmnew <- intron.loss2(gmnew, b, nrt);
	if (i %% 1 == 0 ) {
		h <- statgene(gmnew$reloc)
		hist(h$rel.loc, 100, main='', xlab='', ylab='', xlim=c(0,1), freq=F)
		hh <- hist(h$num.intron, 100, main='', xlab='', ylab='')
		avg.intron <- mean(h$num.intron)
		text(max(hh$breaks)*0.5, max(hh$count)*0.9, 
			labels=paste(i*nrt, "\n", sprintf("%1.2f",avg.intron)), 
			cex=0.9, adj=0)
	}
}
		
		
# loc contains for each gene the number of introns
# as a density function, b is restricted to
# 0 <= b <= 1 if 3' end gets more loss than 5' end
# in another word to have the slope >= 0

inloss.pdf <- function(b,x) {
	2*(1-b)*x+b
}

plot(NULL, xlim=c(0,1), ylim=c(0,2), xlab='Intron Relative Location 5\' to 3\'',
	ylab='')
b <- 0
abline(b, 2*(1-b))
text(0.15,0.2, labels='b=0')
b <- 0.5
abline(b, 2*(1-b), col='red')
text(0.15,0.56, labels='b=0.5')
b <- 1 
abline(b, 2*(1-b), col='green')
text(0.15, 1.1, labels='b=1')

##-- similuation of intron loss from 3' end
## simulated intron loss probility density function
b <- 0
nrt <- 1000
# number of RT event in the whole genome
U <- runif(nrt)
x <- (sqrt(b^2+4*(1-b)*U) - b)/(2*(1-b))
hist(x, 50)

hist(intron.loc, 100)

-- should use sampling with repeats




#-- sample with repeat
b <- 0.9

b <- 0.1
nrt <- 1000
numexon <- NULL
for (i in 1:nrt) {
	pos <- floor(runif(1, 1, ng)+0.5)
	stop.at <- (sqrt(b^2+4*(1-b)*runif(1))-b)/(2*(1-b))
	if (loc[[pos]][1] != 0) {
		tmp <- loc[[pos]]
		if (length(tmp[tmp < stop.at]) < 1) {
			loc[[pos]] <- c(0);
			numexon <- c(numexon, 1);
		}
		else {
			loc[[pos]] <- tmp[tmp<stop.at]
			numexon <- c(numexon, length(loc[[pos]])+1)
		}
	}
}
# now plot loc
N <- length(loc)
h <- numeric()
for (i in 1:N) {
	h <- c(h, loc[[i]])
}
hist(h[h>0], 100, freq=F)
#-- this is very similar to the one observed with Sporo1
hist(numexon, 100)


hist(h, 20)
