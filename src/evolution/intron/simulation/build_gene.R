# function to simulate genes given the distribution of 
# number of introns and the length of exon both obey
# gamma distribution
# I have estimated these parameters based on the 16 fungal genome
# now we have much better EST-annotated genomes, I should
# try to use these genomes.


#ng <- 100
#alpha <- 2.4232631
#lambda <- 0.3130604
#conserve.frac <- 0.2
#conserve.per.gene <- 0.8
#param <- c(ng, alpha, lambda)
#exlen.param.file <- "gamma_param_exlen.tab"
# the table contain up to 70 exons, above it use 70-param
#conserve.param <- c(conserve.frac, conserve.per.gene)

# given a vector of exon lengths, compute relative intron location
# if no intron then return numeric() which
# is an empty vector of length zero!
compute.ril <- function(exons) {
    ril <- numeric(0);
    numintron <- length(exons)-1
    if (numintron > 0) {
        ril <- numeric(numintron);
        gl <- sum(exons); rl <- 0; 
        for (j in 1:numintron) { 
            rl <- rl + exons[j]/gl
            ril[j]=rl
        }
    } 
    return(ril)
}

# 1. Version one is more realistic, because it also simulate
#    The edge effect of exons
#
# Parameters:
#   ng: numer of genes to make
#   nipar: number of introns gamma distribution parameter
#         vector(alpha, lambda)
#   elpf: exlength.parameter.file
#     elpf has the following columns
#     nx    alpha    lambda    mle.alpha    mle.lambda
#     1    1.62903816967233    0.00138031349487282    1.99070381    0.001686759
#     2    1.588107630887    0.00280046638338852    2.162184    0.003812792
#
#   csvp: vector(fraction.of.genes.with.conserved.introns, 
#                fraction.of.introns.conserved.per.gene)
#
### result list(reloc, conserve, exonl)
#     reloc: relative location of intron
#   conserve: conserved intron positions from small to large
#   exonl: exon length
#
# these three properties are parallel for genes and not grouped into 
# individual gene
build.gene <- function(ng, nipar, csvp, elpf, scalel) {
    exlen.param <- read.table(elpf, header=T)
    exlen.param$mle.alpha <- exlen.param$mle.alpha*scalel
    exl <- list(); rel <- list(); csip <- list();
    # number of introns for each gene
    numin <- round(rgamma(ng, nipar[1], nipar[2]))
    for (i in 1:ng) {
        nx.lookup <- numin[i]+1;
        if (numin[i]+1>70) { nx.lookup=70; }
        # simulate one gene
        gp <- exlen.param[nx.lookup,][c(4,5)]
        exons <- floor(rgamma(numin[i]+1, gp[[1]], gp[[2]])+0.5)
        # accumulate result from one gene
        # 1) exon length
        exl <- c(exl, list(exons))
        # 2) relative intron location RIL
        onegene.reloc <- compute.ril(exons);
        rel <- c(rel, list(onegene.reloc)); 
        # 3) conserved intron position
        ci <- numeric()
        if (runif(1) < csvp[1]) { # toss a die
            csvn <- round(csvp[2]*numin[i])
            if (csvn>0) {
                ci <- sort(sample(seq(1,numin[i]), csvn, replace=F))
            }
        }
        csip <- c(csip, list(ci))
    }
    return(list(reloc=rel, conserve=csip, exonl=exl))
}

# 2. Version 2 is simplified, only simulate the intron locations
# whih is uniformation distribution

############################### build genes ########################################
## simulation of intron loss
# it seems that only a smaller fraction of the genes
# are resistatn to gene loss, 80% of the introns
# should be preserved
#-- we assume the Sporolomyces lost 0.5 intron on average
# the number of introns per gene

## the gene is a list containing two objects:
#  1. relative location of each intron
#  2. conserved introns as intron location (1-based index)
# this function will build this object

#alpha <- 2.4232631
#lambda <- 0.3130604
# ng: number of genes to make
# nipar: number of introns gamma distribution parameter
#        vector(alpha,lambda)
# csvp: conserved intron parameter
#         vector(fraction of genes with conserved introns,
#               fraction of introns conserved of a gene)
# this version should not be used in the future.
build.gene.uniform <- function(ng, nipar, csvp) {
# fixed parameters
    numintron <- floor(rgamma(ng, nipar[1], nipar[2])+0.5)
    rel <- list(); consvi <- list();
    for ( i in 1:ng) {
        ci <- 0; v <- 0
        if (numintron[i] > 0) {
            v <- sort(runif(numintron[i]))
            if (runif(1) < csvp[1]) {
                csvn <- floor(csvp[2]*numintron[i]+0.5)
                ci <- sort(sample(seq(1,numintron[i]), csvn, replace=F))
            }
        }
        rel <- c(rel, list(v)); 
        consvi <- c(consvi, list(ci))
    }
    list(reloc=rel, conserve=consvi)
}


# version 3
# this is the same as version 1, but used the 
# fixed parameter for ancient gene features
# both exon length and number of introns are
# from gamma distribution
#
# Parameters:
#   ng: numer of genes to make
##  nipar: number of introns gamma distribution parameter
#         vector(alpha, lambda)
#   elpar: exlength.parameter (alpah, lambda)
#   csvp: vector(fraction.of.genes.with.conserved.introns, 
#                fraction.of.introns.conserved.per.gene)
### result list(reloc, conserve, exonl)
#     reloc: relative location of intron
#   conserve: conserved intron position
#   exonl: exon length
build.fix.gene <- function(ng, nipar, csvp, elpar) {
    #exlen.param <- read.table(elpf, header=T)
    exl <- list(); rel <- list(); csip <- list();
    # number of introns for each gene
    # round by 0.5
    numin <- round(rgamma(ng, nipar[1], nipar[2]))
    for (i in 1:ng) {
        exons <- round(rgamma(numin[i]+1, elpar[1], elpar[2]))
        # accumulate result from one gene
        # 1) exon length
        exl <- c(exl, list(exons))
        # 2) relative intron location RIL
        #    This is derived feature of 1)
        rel <- c(rel, list(compute.ril(exons))); 
        # 3) conserved intron position index
        ci <- 0;
        if (runif(1) < csvp[1]) {
            csvn <- floor(csvp[2]*numin[i]+0.5)
            if (csvn>0) {
                ci <- sort(sample(seq(1, numin[i]), csvn, replace=F))
            }
        }
        csip <- c(csip, list(ci))
    }
    list(reloc=rel, conserve=csip, exonl=exl)
}

# given the genes, it produces statistics for
#   1. relative intron location
#   2. exon length distribution
#   3. number of introns per gene distribution
# a big vector of features for each gene
# all three are required by the paper.
# will not work with genes with 2 feature only
flatten.gene.feat <- function(genes) {
    N <- length(genes$reloc);
    loc <- numeric(); numin <- numeric();
    exl <- numeric();
    for (i in 1:N) {
        EPG <- length(genes$exonl[[i]])
        if (EPG > 1) { # RIL make sense
            loc <- c(loc, genes$reloc[[i]])
            if (length(genes$reloc[[i]])+1 != EPG) {
                print(paste('exons', genes$exonl[[i]]))
                print(paste('ril', genes$reloc[[i]]))
                stop("EPG not matching numin inside flatten.gene.feat")
            }
        }
        numin <- c(numin, EPG-1)
        exl <- c(exl, genes$exonl[[i]])
    }
    return(list(rel.loc=loc, num.intron=numin, exlen=exl))
}

mean.flat.gene <- function(flat.gene) {
    return( c(mean(flat.gene$rel.loc), 
        mean(flat.gene$exlen),
        mean(flat.gene$num.intron)) )
}

# input object flatgenome is the flattened features
# returned from flatten.gene.feat
# also return (mean ril, mean exon length, mean num intron per gene)
# rc is the round count, ancestor has zero
drawhist <- function(gs, rc) {
    mean.feat <- mean.flat.gene(gs)
    #print(paste('mean features of genome:', mean.feat))
    mril <- mean.feat[1]; mexl <- mean.feat[2];
    mipg <- mean.feat[3]
    #par(mfrow=c(3*numr,numc), mar=c(4,6,2,2), mex=0.5, oma=c(3,1,1,1))
    r <- hist(gs$rel.loc, 100, main='', xlab='Relative Intron Location',
        ylab='')
    ps=0.9
    if (rc == 0) { ps=0.5 }
    text(0.85, ps*max(r$count), paste('Avg = ', round(mril,3)), col='red')
    r <- hist(gs$exlen, 500, main='', xlab='Exon Length', ylab='',
        xlim=c(0, 1000))
    #text(0.5*max(gs$exlen), max(r$count), 
    text(500, 0.9*max(r$count), paste('Avg =', round(mexl, 1)), col='red')
    r <- hist(gs$num.intron, 100, main='', xlab='Number of Introns per Gene', 
        ylab='')
    text(0.5*max(gs$num.intron), 0.9*max(r$count), 
        paste('Avg =', round(mipg, 1), ' round =', rc), col='red')
    return(c(mril, mexl, mipg))
}

# directly from genome
# input genome drawn on three pannels on the same figure
drawhistGenome <- function(genome) {
    gs <- flatten.gene.feat(genome)
    mean.feat <- mean.flat.gene(gs)
    mril <- mean.feat[1]; mexl <- mean.feat[2];
    mipg <- mean.feat[3]
    par(mfrow=c(3,1), mar=c(4,6,1,2), oma=c(3,2,2,1), mex=0.8, cex=0.8)
    r <- hist(gs$rel.loc, 100, main='', xlab='Relative Intron Location',
        ylab='')
    text(0.5, 0.5*max(r$count), cex=1.1, paste('Avg = ', round(mril,3)), col='red')
    r <- hist(gs$exlen, 500, main='', xlab='Exon Length', ylab='',
        xlim=c(0, 1000))
    #text(0.5*max(gs$exlen), max(r$count), 
    text(500, 0.9*max(r$count), paste('Avg =', round(mexl, 1)), col='red')
    r <- hist(gs$num.intron, 100, main='', xlab='Number of Introns per Gene', 
        ylab='')
    text(0.5*max(gs$num.intron), 0.9*max(r$count), 
        paste('Avg =', round(mipg, 1)), col='red')
    return(c(mril, mexl, mipg))
}

################# testing the two functions ########################
createFungalAncestor <- function(ng, alpha, lambda, scalel, conserve.param) {
    #conserve.frac <- 0.2
    #conserve.per.gene <- 0.8
    ### combined parameters
    numintron.par <- c(alpha,lambda)
    exlen.param.file <- "gamma_param_exlen.tab"
    #conserve.param <- c(conserve.frac, conserve.per.gene)

    gmodel <- build.gene(ng, numintron.par, conserve.param, 
        exlen.param.file, scalel);
    #simple.gmodel <- build.gene.uniform(ng, numintron.par, conserve.param);
    #gs <- gene.stat(gmodel);
    #drawhist(gs)
    return(gmodel)
}

createEukaryoteAncestor <- function(ng, nipar, exlpar, consvpar) {
    gmodel <- build.fix.gene(ng, nipar, consvpar, exlpar);
    #gs <- flatten.gene.feat(gmodel);
    #drawhist(gs)
    return(gmodel)
}



#ng <- 10000
# 1. Roughly value obtained by RIL method
#alpha <- 2.4232631
#lambda <- 0.3130604
#conservepar <- c(0.2, 0.8)
#fmodern.genome <- testFungalAncestor(ng, alpha, lambda, 1, conservepar)
# simulation result:
# RIL: 0.499, Exlen: 213.7 nt, IPG: 7.666

# 2. By RT method about 10 EPG
# since exon length will be to too long
# I will scale the alpha by 0.65
#alpha <- 2.297
#lambda <- 0.255
#fancient.genome <- testFungalAncestor(ng, alpha, lambda, 0.655, conservepar)
# simulation result:
# RIL: 0.500, Exlen: 133.4 nt, IPG: 9.0
# CDS length 1200


# 3. Ancestor of all eukaryotes

# ancestral gene
#exlen.param <- c(1.1408, 0.01755)
# mean= 65.00285 nt exon length
#numin.param <- c(4.003, 0.223)
# mean = 17.95 exons per gene
#consv.param <- c(0.2, 0.8)
#eancient.genome <- testEukaryoteAncestor(ng, numin.param, 
#    exlen.param, consv.param)

# after the run this program produces three genomes
