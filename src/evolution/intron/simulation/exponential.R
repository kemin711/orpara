# similuation with linear over the range [0,1] 
# is not realistic, so I am using the result
# from RT enzymology studies
# it was found that RT has dissociation rate
# 0.0026/nt and for retroviruses it is in
# the range of 0.0054/nt. retrovirueses have
# higher dissociation rates
# now we also need to model the length of the 
# gene models.

-- for Sporo1 whose parameter is simulated by gamma
-- distribution: alpha=1.9232631, lambda <- 0.3130604
-- This is the number of exons distribution
-- we assume this as the ancestral state
-- assume 10,000 genes in the genome

alpha <- 1.9232631
lambda <- 0.3130604
x <- c(seq(0,10, 0.01), seq(10.2, 70, 0.5))
plot(x, dgamma(x, alpha, lambda))
lines(x, dgamma(x, alpha+0.5, lambda))

ng <- 10000
#-- we assume the Sporolomyces lost 0.5 intron on average
# number of introns per gene
numin <- floor(rgamma(ng, alpha+0.5, lambda)+0.5)
# then we simulate the exons with statistical features
# from all fungal genomes. Exon length shows a gamma
# distribution.  With average exon length v.s. number
# of exons per gene assuming a complicated relationship
exlen.param <- read.table("gamma_param_exlen.tab", header=T)

###### RT Process #########
#poff <- 0.0026
poff <- 0.0021


# plot the intron location distribution
# gg is the all genes: list of list of vector
relocation <- function(gg) {
	loc <- numeric()
	for (i in 1:length(gg)) {
		g <- gg[[i]]
		if (length(g)>1) {
			total <- 0
			for (j in 1:(length(g)-1)) {
				total <- total + g[j]/sum(g)
				loc <- c(loc, total)
			}
		}
	}
	loc
}

### intron loss function, g is a gene with vector structure
# g (exon1.len, exon2.len, ...)
# trlen is the length of the RT generated transcript
intron.loss <- function(g, trlen) {
	i <- length(g);
	if (i == 1) {
		new.g <- g;
		return(new.g)
	}
	exend.sum <- g[i];
	#print(paste('initial sum', exend.sum))
	while (exend.sum < trlen & i > 1) {
		i <- i - 1;
		exend.sum <- exend.sum + g[i];
		#print(paste('sum:', exend.sum, ' index:', i))
	}
	if (i>1) {
		new.g <- c(g[seq(1,i-1)], exend.sum)
	}
	else { new.g <- exend.sum }
	#print(g)
	#print(paste('trlen', trlen))
	#print(new.g)
	return(new.g)
}

# get information about exon length of all the 
# genes
exon.length <- function(gg) {
	exlen <- numeric()
	for (i in 1:length(gg)) {
		exlen <- c(exlen, gg[i][[1]])
	}
	exlen
}

number.exon <- function(gg) {
	nx <- numeric()
	for (i in 1:length(gg)) {
		nx <- c(nx, length(gg[i][[1]]))
	}
	nx
}

### 1. build genes #####
genes <- list()

for (i in 1:length(numin)) {
#for (i in 1:100) { # for testing
	while (numin[i] > 70) {  # make sure less than 70
		numin[i] <- floor(numin[i]*runif(1))
	}
	gp <- exlen.param[numin[i]+1,][c(4,5)]
	exons <- list( floor(rgamma(numin[i]+1, gp[[1]], gp[[2]])+0.5) )
	genes <- c(genes, exons)
}
reloc <- relocation(genes)
par(mfrow=c(4,4))
hist(reloc, 100, freq=F)

########################################
### build gene for linear loss function 
## in linear_conserved.R intron.loss2()




### 2. we assume that 100 RT events happened to the genome


num.rt <- 100
simcount <- 100

tmp <- genes
for (s in 1:simcount) {
	pick <- floor(runif(num.rt, 1, ng)+0.5)
	rtlen <- rexp(num.rt, poff)
	for (p in 1:num.rt) {
		new.exons <- intron.loss(tmp[pick[p]][[1]], rtlen[p])
		tmp[pick[p]][[1]] <- new.exons
	}
}
# for plotting use relocation	
reloc <- relocation(tmp)
hist(reloc, 100, freq=F)	

el <- exon.length(tmp)
hist(el, 100)

numex <- number.exon(tmp)
hist(numex, 100)

# the result of the simulation is not exact as 
# expected, there are too many single exon genes
# after intron loss step.
# need to add 3'-polyA factor, 5'-nuclease digestion
# length of the DNA and the recombination opportunity
# probility


