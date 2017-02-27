enl <- read.table("exnum_avgexlen.tab", header=T)
names(enl)
[1] "db"  "proteinid" "numex"     "avgexlen"

plot(enl$numex, enl$avgexlen)

mean(enl$avgexlen)
# 529
sd(enl$avgexlen) # 621
median 322

enl.noextreme <- subset(enl, avgexlen<10000)     
mean(enl.noextreme$avgexlen)
# 527, not much difference

par(mfrow=c(3,3))
xmax=c(8000,4000,3000,2000,2000,1500,1500, 1500,1000,
	1000,900,900, 900,800,800, 800,800, 800);

xrange <- range(enl.noextreme$avgexlen)
x <- seq(0,xrange[2], 0.5)
for (i in 9:17) {
for (i in 1:9) {
	tmp <- subset(enl.noextreme, numex==i & avgexlen < xmax[i])
	avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
	count <- length(tmp$avgexlen)
	lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
	gamma.val <- dgamma(x, alpha, lambda)
	h <- hist(tmp$avgexlen, 100, main=i)
	top <- max(h$counts)
	fac <- max(h$counts)/max(gamma.val)
	lines(x, 0.8*gamma.val*fac, col='red')
	text(xmax[i]-200, top*0.9, labels=paste(count, "\n", alpha, "\n", lambda), adj=1)
}

estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 1.99070381)
alpha <- 1.99070381
lambda <- alpha/avg.exl

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, gamma.val*fac, col='red')

# area should match to max
# estimate the formular for mean exon length against
# the number of exons (or introns) per gene
# the raw data is taken from exnumlen_summary.tab
enls <- read.table("exnumlen_summary.tab", header=T);

names(enls)
[1] "numex"   "meanlen" "sdlen"   "count"

# searching for the best parameters
miny <- 1000000
arange <- seq(900, 1400, 10)
brange <- seq(-0.8, -0.2, 0.01)
crange <- seq(-2, -1.1, 0.01)
drange <- seq(160, 240, 4)
[1] 1060.00   -0.53   -1.42  204.00
[1] 1055.00   -0.53   -1.42  204.00
[1] 1056.00   -0.53   -1.42  204.00
[1] 1056.00   -0.53   -1.42  204.00
[1] 1056.000   -0.530   -1.419  204.000
[1] 1056.000   -0.530   -1.419  204.000
[1] 1057.0000   -0.5304   -1.4192  204.0000

-- use log distince get more weight for smaller 
# numbers
miny <- 1000000
arange <- seq(1056, 1059, 1)
brange <- seq(-0.535, -0.525, 0.0001)
crange <- seq(-1.42, -1.411, 0.0001)
drange <- seq(203, 205, 1)
for (a in arange) {
	for (b in brange) {
		for (c in crange) {
			for (d in drange) {
eqn <- function(x) {
	a*exp(b*x)+c*x+d
}
				res <- sum(abs(log(eqn(enls$numex))-log(enls$meanlen)))
				if (res < miny) {
					miny <- res
					minp <- c(a, b, c, d)
				}
			}
		}
	}
}
print(minp)
print(miny)
a <- minp[1]; b <- minp[2]; c <- minp[3]; d <- minp[4]
eqn <- function(x) {
	a*exp(b*x)+c*x+d
}
################## linear penalty ############
972.00  -0.79  -1.53 208.00
975.00  -0.78  -1.47 205.00
975.00  -0.78  -1.47 205.00
975.000  -0.784  -1.463 205.000

# use intron number = exnum - 1
miny <- 100000
minp <- NULL
arange <- seq(970, 980, 1)
brange <- seq(-0.79, -0.77, 0.001)
crange <- seq(-1.5, -1.45, 0.001)
drange <- seq(203, 207, 1)
for (a in arange) {
	for (b in brange) {
		for (c in crange) {
			for (d in drange) {
eqn <- function(x) {
	a*exp(b*x)+c*x+d
}
				res <- sum(abs(eqn(enls$numex-1)-enls$meanlen))
				if (res < miny) {
					miny <- res
					minp <- c(a, b, c, d)
				}
			}
		}
	}
}
print(minp)
[1] 975.000  -0.784  -1.463 205.000
minp <- c(975.000, -0.784, -1.463,205.000)
print(miny)
exlnfunc <- function(p,x) {
	p[1]*exp(p[2]*x)+p[3]*x+p[4]
}
numin <- enls$numex-1
par(mfrow=c(2,1), mex=0.5, mar=c(4,5,1,1), oma=c(1,1,1,1), cex=0.95)
plot(numin, enls$meanlen, xlab='', ylab='Average Exon Length')
#lines(enls$numex, exlnfunc(minp,enls$numex))
lines(numin, exlnfunc(minp,numin))
text(30, 1000, "975*exp(-0.784*x) - 1.463*x + 205")
plot(numin, (numin+1)*enls$meanlen, xlab="Number of Introns",
	ylab='Gene Length')
lines(numin, (numin+1)*exlnfunc(minp,numin))
 
###### Final Results with log penalty  #################
a <- 1057; b <- -0.5304;  c <- -1.4192; d <- 204

## where x is number of exons per gene
avg.exlen <- 1057*exp(-0.5304*x) -1.4192*x + 204
eqn <- function(x) {
	a*exp(b*x)+c*x+d
}
plot(enls$numex, enls$meanlen)
x <- enls$numex
lines(x, eqn(x))

