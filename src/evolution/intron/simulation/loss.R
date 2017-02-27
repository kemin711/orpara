#-- for Sporo1 whose parameter is simulated by gamma
#-- distribution: alpha=1.9232631, lambda <- 0.3130604
#-- This is the number of exons distribution
#-- we assume this as the ancestral state
#-- assume 10,000 genes in the genome

# number of introns per gene
# the parameters for Sporo1, which has on average, 6.2 introns
# The conserved gens of Lacbi1 has 6.9 introns, so we we assume
# the ancestor has 7.8 introns
#alpha <- 1.9232631
#lambda <- 0.3130604
#x <- c(seq(0,1, 0.02), seq(1.02, 70, 0.5))
#plot(x, dgamma(x, alpha, lambda))
#lines(x, dgamma(x, alpha+0.5, lambda))

##       1          3          5         7           9        11   12
#tf <- c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE);
# number of F before index P
numberf <- function(TF, P) {
    tmpx <- TF[1:P-1]
    length(tmpx[tmpx == F])
}

# tf will be modified and returned
# parameter
#  tf: saved intron T/F vector. Usually 5' to RT stop position
#  oldpos: conserved introns position, integer
newconsvpos <- function(tf, oldpos) {
    newpos <- numeric()
    for (i in 1:length(oldpos)) {
        tf[oldpos[i]] <- T; # if preserved always preserved
        # shift index by the number of introns removed
        # 5' to this location
        newpos <- c(newpos, oldpos[i] - numberf(tf,oldpos[i]))
    }
    list(tf,newpos)
}

# return n numbers between 0 and 1
# with density f(x)=ax+b or f(x)=2(1-b)x+b
# 0<= b <= 1
dintronloss <- function(b, n) {
    U <- runif(n)
    X <- 0.5*(sqrt(b^2+4*(1-b)*U)-b)/(1-b)
    return(X)
}


### testing
#tf <- c(    TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE);
#exons <- c(10, 20,   30,  40,   50,  60,  70,  80,   90,  100, 110, 120, 130);
## count how many FALSE in the vector
# not sure this is needed because it is only one line
#loss.count <- function(tf) {
#    length(tf[tf == FALSE])
#}

## givenn a list of exons, and save/remove boolean vector
# parameters
#   exlen: exon length vector
#   psvtf: preserve T/F
# return a new exon length vector after intron loss
merge.exon <- function(exlen, psvtf) {
    if (length(exlen) != length(psvtf)+1) {
        print(exlen)
        print(psvtf)
        stop("exon length vector should be 1-longer than intron preserve vector");
    }
    newex <- numeric();
    i=1; ni <- length(psvtf);
    while (i <= ni) {
        nwx <- exlen[i];
        #print(paste('ni=', ni, 'i:', i, 'psvtf[i]:', psvtf[i]));
        while (i <= ni & psvtf[i] == FALSE) {
            nwx <- nwx + exlen[i+1]; i <- i+1;
        }
        i <- i+1
        newex <- c(newex, nwx)
    }
    if (i == ni+1) {
        newex <- c(newex,exlen[i]);
    }
    #print('old exons: '); print(exlen)
    #print('after merge: '); print(newex)
    return(newex);
}

# testing
#mgd <- merge.exon(exons,tf)
#print(mgd)

### new  version that will use gmodel as input
# it preserved conserved introns 80% of the introns
# the same RT could hit the gene twice or more!
# we sample with replacement.
# exon length seems to be not updated
lose.intron <- function(gm, b, nrt) {
    numgene <- length(gm$reloc); 
    #hasexon <- length(gm)>2;
    # same gene could be affected twice
    gi <- sample(seq(1:numgene), nrt, replace=TRUE)
    X <- dintronloss(b, nrt)
    for (i in 1:nrt) {
        pos <- gi[i]; 
        #print(paste('RT happend to gene ', pos))
        EPG <- length(gm$exonl[[pos]])
        if (EPG > 1) {  # has intron, so test lose
            savedidx <- gm$reloc[[pos]] < X[i] # boolean vector
            # check any one is conserved
            if (length(gm$conserve[[pos]]) > 0) { 
                tmpl <- newconsvpos(savedidx, gm$conserve[[pos]]);
                gm$conserve[[pos]] <- tmpl[[2]]; # new conserved position
                savedidx <- tmpl[[1]]; # Logical Vector
            }
            #gm$reloc[[pos]] <- tmp[savedidx] # this does not work
            # res is the new relative location
            res <- gm$reloc[[pos]][savedidx];
            if (length(res) == 0) { gm$reloc[[pos]] <- list(0); }
            else { gm$reloc[[pos]] <- res; }
            if (length(savedidx[savedidx==FALSE])>0) {
               # applied to gene models with three features
                #print(paste("trying to merge exons", pos))
                gm$exonl[[pos]] <- merge.exon(gm$exonl[[pos]], savedidx)
            }
        }
    }
    return(gm)
}

# N: number of rounds of simulation
# nrt: number of RT events in on round of simulation
# gm: input genome with a pool of genes
# b: RT loss linear density function parameter [0,1]
#
loss <- function(gm, b, nrt, N) {
    par(mfrow=c(8,3), mar=c(2,2,2,2), mex=0.5)
    gm.stat <- flatten.gene.feat(gm)
    mfeat <- drawhist(gm.stat, 0)
    #mfeat <- mean.flat.gene(gm.stat)
    #par(mfrow=c(3,1), mar=c(2,2,4,2), mex=0.9)
    res <- data.frame(round=0, numrt=0,
            ril=mfeat[1], exonlen=mfeat[2], IPG=mfeat[3]);
    print(paste('b=',b))
    ni2loc <- numeric()
    total.rt <- 0
    drawat <- c(1, 2, 4, 8, 16, 32, 64)
    
    for (i in 1:N) {
        #print(paste('lost round: ', i))
        gm <- lose.intron(gm, b, nrt);
        gm.stat <- flatten.gene.feat(gm)
        mfeat <- mean.flat.gene(gm.stat)
        total.rt <- total.rt + nrt
        tmpres <- data.frame(round=i, numrt=total.rt,
            ril=mfeat[1], exonlen=mfeat[2], IPG=mfeat[3]);
        res <- rbind(res, tmpres)
        if ( any(i == drawat)) {
		    print(paste('lost round: ', i))
            drawhist(gm.stat, i)
        }
    }
    print(res)
    write.table(res, file='losstrack.tab', row.names=F,
        quote=F, sep="\t")
    return(res)
}

lossNoDraw <- function(gm, b, nrt, N) {
    #par(mfrow=c(8,3), mar=c(2,2,2,2), mex=0.5)
    gm.stat <- flatten.gene.feat(gm)
    #mfeat <- drawhist(gm.stat, 0)
    mfeat <- mean.flat.gene(gm.stat)
    #par(mfrow=c(3,1), mar=c(2,2,4,2), mex=0.9)
    res <- data.frame(round=0, numrt=0,
            ril=mfeat[1], exonlen=mfeat[2], IPG=mfeat[3]);
    #print(paste('b=',b))
    ni2loc <- numeric()
    total.rt <- 0
    #drawat <- c(1, 2, 4, 8, 16, 32, 64)
    
    for (i in 1:N) {
        #print(paste('lost round: ', i))
        gm <- lose.intron(gm, b, nrt);
        gm.stat <- flatten.gene.feat(gm)
        mfeat <- mean.flat.gene(gm.stat)
        total.rt <- total.rt + nrt
        tmpres <- data.frame(round=i, numrt=total.rt,
            ril=mfeat[1], exonlen=mfeat[2], IPG=mfeat[3]);
        res <- rbind(res, tmpres)
        #if ( any(i == drawat)) {
		  #  print(paste('lost round: ', i))
        #    drawhist(gm.stat, i)
        #}
    }
    #print(res)
    #write.table(res, file='losstrack.tab', row.names=F,
    #    quote=F, sep="\t")
    return(res)
}

