#Simulate the RT fall off probility function
lambda=0.0026
par(mfrow=c(1,1), mar=c(5,5,3,1))
curve(dexp(x,lambda), 0, 2000, 1000, xlab='Nucleotides', ylab='Density')
curve(dexp(x,0.001), 0,2000, add=T, col='green')
curve(dexp(x,0.0005), 0,2000, add=T, col='red')
xloc=1600; yloc=0.0023;
text(xloc-5, yloc+0.00012, 'lambda')
text(xloc,yloc, '0.0026', adj=0);
text(xloc,yloc-0.00012, '0.001', col='green', adj=0);
text(xloc,yloc-0.00024, '0.0005', col='red', adj=0);

text(1500,0.0020, paste('Black: 0.0026', "\n", 'Green:0.001', "\n",
	'Red:0.0005'), adj=0.2);

