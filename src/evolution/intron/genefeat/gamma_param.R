## estimate the gamma parameters for exon length against different
## number of exons, each is a gamma distribution.
# I have estimated it manually for exon number
# from 1 to 70. For larger number of exons we usually don't have
# enought number of exon, especially for MLE metheds where
# it performed poorly after 25, we used the Moment method.

enl <- read.table("exnum_avgexlen.tab", header=T)
enl.noextreme <- subset(enl, avgexlen<10000)
maxlen <- 10000

### method of Moments ####################
lambda = avg/var
alpha = avg^2/var

names(enl.noextreme)
[1] "db"        "proteinid" "numex"     "avgexlen"




range(enl.noextreme$numex)
res <- numeric()

for (i in 1:max(enl.noextreme$numex)) {
	print(i)
	tmp <- subset(enl.noextreme, numex==i)
	if (length(tmp$avgexlen) > 0) {
		xl <- tmp$avgexlen
		avg.exl <- mean(xl); var.exl <- var(xl)
		lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
		res <- c(res, i, alpha, lambda)
	}
	else { print("empty") }
}

param.moment <- matrix(res, ncol=3, byrow=T)

write.table(param.moment, 'gamma_param_exlen.tab', quote=F, 
	sep="\t")


par(mfrow=c(1,1))
i <- 1
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0, 2400))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 1.99070381)
alpha <- 1.99070381
lambda <- alpha/avg.exl
lambda
# 0.00168676

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

#########
i <- 2
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0, 5400))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 2.162184)
alpha <- 2.162184 
lambda <- alpha/avg.exl
lambda
0.003812792

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

################### 3 ##################
i <- 3
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0, 6600))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 2.26468426)
alpha <- 2.26468426 
lambda <- alpha/avg.exl
lambda

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 4
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0, 3000))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 2.47752907)
alpha <- 2.47752907 
lambda <- alpha/avg.exl
lambda
# 0.007717834

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 5
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0, 1100))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 2.8483275)
alpha <- 2.8483275 
lambda <- alpha/avg.exl
lambda

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.91*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 6
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0, 2000))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 2.9668215)
alpha <- 2.9668215 
lambda <- alpha/avg.exl
lambda

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.85*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 7
tmp <- subset(enl.noextreme, numex==i & avgexlen < maxlen)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 3.3035991)
alpha <- 3.3035991
lambda <- alpha/avg.exl
lambda

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 8
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 3.551456)
alpha <- 3.551456
lambda <- alpha/avg.exl
lambda


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 9
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(1000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
3.85038752243195   0.0187367749038938


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 10
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 3.99915)
alpha <- 3.99915
lambda <- alpha/avg.exl
lambda

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

par(mfrow=c(1,1))
i <- 11
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(600, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 4.224548)
alpha <- 4.224548
lambda <- alpha/avg.exl
lambda
# 0.02228848

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 12
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(600, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 4.22305)
alpha <- 4.22305
lambda <- alpha/avg.exl
lambda

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.70*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 13
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(600, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 4.42676)
alpha <- 4.42676
lambda <- alpha/avg.exl
lambda
# 0.02380656

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 14
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(600, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 5.07448)
alpha <- 5.07448
lambda <- alpha/avg.exl
lambda
# 0.0290384

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.70*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)


i <- 15
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 3.7522)
alpha <- 3.75223
lambda <- alpha/avg.exl
lambda
# 0.02065346


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 16
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(500, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 4.18268)
alpha <- 4.18268
lambda <- alpha/avg.exl
lambda
# 0.02388844


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 17
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 4.18679)
alpha <- 4.18679
lambda <- alpha/avg.exl
lambda


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 18
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
error <- 0.00001
val <- 10
step <- 0.1
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
}
lambda <- alpha/avg.exl
print(paste(alpha, "\t", lambda))
4.33980989276027   0.0253612800645335


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 19
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, " ", lambda))
4.17071774166635 0.0224068909905155


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 20
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 4.66119)
alpha <- 4.66119
lambda <- alpha/avg.exl
lambda


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 21
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
3.89466731174051 0.0216142660888801

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 22
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
3.76927166357573 0.0221413455327209

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 23
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 4000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
5.01401275106742 0.0307460845421494

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 24
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 4000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
5.10779500937041 0.0309417533575529

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 25
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 5.1545)
alpha <- 5.1545
lambda <- alpha/avg.exl
lambda


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 26
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 4000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
5.2162306001017 0.0317142920051124

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 27
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 4000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
5.33522439561174 0.0288875414181187

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(400, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 28
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 4000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
4.15986475712478 0.0247457662075153

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 29
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex==i & avgexlen < 4000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
4.64713279116989 0.0240542740770288

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

### not enough data for MLE
i <- 30
tmp <- subset(enl.noextreme, numex==i & avgexlen < 5000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i, ylim=c(0,3))
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(50, top, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

### canno estimate any more, not enough data!
### after the estimate
estalpha(count, tmp$avgexlen, alpha)
estalpha(count, tmp$avgexlen, 3.865)
alpha <- 3.865
lambda <- alpha/avg.exl
lambda


h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.80*gamma.val*fac, col='red')
text(5000, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

i <- 31-35
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex >= 31 & numex <= 35 & avgexlen < 3000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
5.00278589698288 0.0317955157770044

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)


i <- 36-70
par(mfrow=c(1,1))
tmp <- subset(enl.noextreme, numex >= 36 & avgexlen < 3000)
avg.exl <- mean(tmp$avgexlen); var.exl <- var(tmp$avgexlen);
count <- length(tmp$avgexlen)
lambda <- avg.exl/var.exl; alpha <- avg.exl^2/var.exl;
gamma.val <- dgamma(x, alpha, lambda)
h <- hist(tmp$avgexlen, 100, main=i)
top <- max(h$counts)
fac <- max(h$counts)/max(gamma.val)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.6f",alpha), "\n", 
	sprintf("%1.6f",lambda)), adj=1)
print(paste(alpha, lambda))
3.27274263513026 0.0229888622234749

error <- 0.00001; val <- 10; step <- 0.1;
while (abs(val) > error) {
	while (val > 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha + step;
	}
	step <- step/5
	while (val < 0 & abs(val) > error) {
		val <- estalpha(count, tmp$avgexlen, alpha)
		alpha <- alpha - step;
	}
	step <- step/5
}
lambda <- alpha/avg.exl
print(paste(alpha, lambda))
2.05930263513026 0.0144652451577709

h <- hist(tmp$avgexlen, 100, main=i)
gamma.val <- dgamma(x, alpha, lambda)
lines(x, 0.95*gamma.val*fac, col='red')
text(300, top*0.9, labels=paste(count, "\n", sprintf("%1.8f",alpha), "\n", 
	sprintf("%1.8f",lambda)), adj=1)

# the estimate by MLE is very bad, Moments is better
