enl <- read.table("exnum_avgexlen.tab", header=T)
enl.noextreme <- subset(enl, avgexlen<10000)
gamma.param <- read.table("gamma_param_exlen.tab", header=T)
names(enl.noextreme)
[1] "db"        "proteinid" "numex"     "avgexlen"

## the whole database, show the first 16
intron.numlen <- data.frame(numintron=enl.noextreme$numex-1, avgexlen=enl.noextreme$avgexlen)
par(mfrow=c(4,4), mex=0.6, mar=c(3,3,1,1))
N <- 15
xend <- 10000
grow <- c(1,1,1,1, 1,1,1,1, 1,0.8,0.8,0.8, 0.87,0.75,0.8,0.6)
x <- c(seq(0,500,0.05), seq(500.1, xend, 1))
for (i in 0:N) {
	tmp <- subset(intron.numlen, numintron==i & avgexlen<xend);
	gp <- subset(gamma.param, nx==i+1)
	if (max(tmp$avgexlen) < xend) {
		xend <- max(tmp$avgexlen)
	}
	h <- hist(tmp$avgexlen, 50, xlim=c(0,xend), main='',
		xlab='', ylab='')
	tx <- xend*0.8; ty <- max(h$counts)
    text(tx, ty*0.8, labels=paste("NI:", i, "\n", sprintf("%1.3f", gp$mle.alpha),
		"\n", sprintf("%1.5f", gp$mle.lambda)), adj=1)
	gval <- dgamma(x, gp$mle.alpha, gp$mle.lambda)
	fc <- ty/max(gval)
	lines(x, grow[i+1]*fc*gval, col='blue')
}

