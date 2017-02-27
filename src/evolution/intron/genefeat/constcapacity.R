gnz <- read.table("gensz_numgene.tab", header=T)

names(gnz);

xleft <- min(gnz$genome_size)*0.83
xright <- max(gnz$genome_size)
str <- c("Ng=2.573e-04G + 1.278e+03", "p-val=2.765e-07")

par(mex=0.7)
plot(gnz$genome_size, gnz$numgene, xlab='Genome Size', 
	ylab="Number of Genes", xlim=c(xleft,xright));
for (i in 1:length(gnz$genome_size)) {
	tp <- 2;
	if (gnz[i,2] == 'Sporo1' | gnz[i,2] == 'ustma1' | gnz[i,2] == 'Mycgr1') {
		tp <- 4;
	}
	else if (gnz[i,2] == 'Aspni1') { tp <- 1; }
	else if (gnz[i,2] == 'Trive1') { tp <- 3; }
	text(gnz[i,3], gnz[i,4], gnz[i,2], cex=0.6, pos=tp)
}

seldb <- subset(gnz, dbid != 6 & dbid != 12);

lm.cap <- lm(seldb$numgene ~ seldb$genome_size)
abline(lm.cap, col='red')
text(1.5e+07, 16000, str[1], cex=0.8, pos=4)
text(1.5e+07, 15400, str[2], cex=0.8, pos=4)

