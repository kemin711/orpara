nxl <- read.table("exlen_vs_numex.tab", header=T)
# in this file I grouped exon length against the number of exons
# no grouping was done for the length

numex.max <- max(nxl$numex);

#for (i in 1:numex.max) {

par(mfrow=c(4,4))
for (i in 1:16) {
	ds <- subset(nxl, numex==i);
	plot(ds$exlen, ds$count)
}

## grouped by 10
nxl <- read.table("exlen_vs_numex10.tab", header=T)

numex.max <- max(nxl$numex);

#for (i in 1:numex.max) {

par(mfrow=c(4,4), mex=0.6)
for (i in 1:16) {
	ds <- subset(nxl, numex==i);
	plot(ds$exlen, ds$count, xlab='')
}


## grouped by 20
nxl <- read.table("exlen_vs_numex20.tab", header=T)
sta <- read.table("numex_statexlen.tab", header=T)


numex.max <- max(nxl$numex);

#for (i in 1:numex.max) {

par(mfrow=c(4,4), mar=c(3,3,0,2), mex=0.6, oma=c(2,2,2,2))
for (i in 1:16) {
	ds <- subset(nxl, numex==i);
	plot(ds$exlen*20, ds$count, xlab='', ylab='')
	text(max(ds$exlen*20)*0.7, max(ds$count)*0.8, 
		paste("Ni: ", i, "\n",
		floor(sta[i, 2]+0.5), " +/- ", floor(sta[i,3]+0.5))
		)
}

# obtained by the best parameter
ELN <- function(x) {
	1060.1*exp(-0.7812*x) + -1.8961*x + 198.5;
}
# some variation
ELN <- function(x) {
	1060.1*exp(-0.7822*x) + -1.4961*x + 198.5;
}

par(mfrow=c(2,1), mar=c(2,5,1,2), mex=0.5, oma=c(4,1,1,1))
x <- seq(0,70,0.1)
plot(sta$numex-1, sta$avgexlen, ylab='Exon Length', xlab='')
lines(x, ELN(x))
plot(sta$numex-1, sta$avgexlen*sta$numex, ylab='Gene Length', 
	xlab='Number of Introns')
lines(x, ELN(x)*(x+1))


## estimate the parameters of the equation
miny <- 1000000
minp <- NULL

arange <- seq(500, 1500, 50)
brange <- seq(-2, -0.01, 0.1)
crange <- seq(-3, -0.5, 0.5)
drange <- seq(100, 300, 10)

# round one result:
# 980.00  -0.77  -1.45 207.00
miny <- 100000
minp <- NULL
arange <- seq(600, 1200, 20)
brange <- seq(-1.2, -0.2, 0.1)
crange <- seq(-2.8, -0.8, 0.4)
drange <- seq(120, 280, 10)
#1060.0   -0.8   -2.0  200.0
# 1060.0   -0.8   -1.9  200.0
# 1060.00   -0.75   -1.80  196.00
# 1059.00   -0.78   -1.90  199.00
# 1060.00   -0.78   -1.86  198.00
# 1060.50   -0.78   -1.86  198.00
# 1060.50   -0.78   -1.86  198.00
# 1060.600   -0.780   -1.868  198.000
# 1060.200   -0.781   -1.880  198.400
# 1060.200   -0.781   -1.890  198.400 @ 918
# 1060.2000   -0.7810   -1.8905  198.4000
# 1060.1000   -0.7812   -1.8961  198.5000

miny <- 100000
minp <- NULL
arange <- seq(1060.1, 1060.3, 0.1)
brange <- seq(-0.782, -0.780, 0.0001)
crange <- seq(-1.90, -1.88, 0.0001)
drange <- seq(198.1, 198.6, 0.1)

for (a in arange) {
	for (b in brange) {
		for (c in crange) {
			for (d in drange) {
eqn <- function(x) {
	a*exp(b*x)+c*x+d
}
				res <- sum(abs(eqn(sta$numex-1)- sta$avgexlen))
				if (res < miny) {
					miny <- res
					minp <- c(a, b, c, d);
				}
			}
		}
	}
}
minp
miny

ELN <- function(x) {
	minp[1]*exp(minp[2]*x) + minp[3]*x + minp[4]
}

ELN <- function(x) {
	1060.1*exp(-0.7812*x) + -1.8961*x + 198.5;
}

par(mfrow=c(1,1))
plot(sta$numex-1, sta$avgexlen*sta$numex)
a <- 200
b <- 200
exl <- function(x) {
	a*x + b
}
lines(x, minp[1]*x + minp[2])

miny <- 5000000;
arange <- seq(1,300,10)
brange <- seq(1,300,10)
for (a in arange) {
	 for (b in brange) {
exl <- function(x) { 
	a*x + b;
}
	 	res <- sum(abs(exl(sta$numex-1) - sta$avgexlen*sta$numex))
		if (res < miny) {
			miny <- res
			minp <- c(a,b)
		}
	}
}
miny
minp


par(mfrow=c(1,1))
plot(sta$numex-1, sta$avgexlen*sta$numex)
lines(x, minp[1]*x + minp[2])

#-- all exon length
exl <- read.table("allcdsexlen.tab", header=TRUE)

par(mfrow=c(1,1))
plot(exl$exlen, exl$count)

exl1 <- subset(exl, exlen<500)
plot(exl1$exlen, exl1$count)

exlp0 <- subset(exl1, exlen %% 3 == 0)
exlp1 <- subset(exl1, exlen %% 3 == 1)
exlp2 <- subset(exl1, exlen %% 3 == 2)

plot(NULL, xlim=c(0,500), ylim=c(0,3000), ylab='Count', xlab='Exon Length')
points(exlp0$exlen, exlp0$count, pch=1, col='red');
points(exlp1$exlen, exlp1$count, pch=2, col='green');
points(exlp2$exlen, exlp2$count, pch=3, col='blue');
legend(400,2700, c("0  280067", "1  206910", "2  206494"), pch=c(1,2,3),
	col=c("red", "green", "blue"))

text(300, 2500, "Phase 0 280067\n 1 206910\n 2 206494", adj=1)


