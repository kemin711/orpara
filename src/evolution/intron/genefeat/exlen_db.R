/* exlength distribution grouped by each genome 
*/

eld <- read.table("cdsexlen_db.tab", header=TRUE)
lenlimit <- 600;
dbs <- levels(eld$db);
par(mfrow=c(4,4), mex=0.7, mar=c(2,3,2,1))
for (d in dbs) {
	dta <- subset(eld, db==d & exlen < lenlimit);
	plot(dta$exlen, dta$count, xlab='', ylab='');
	text(400, max(dta$count)*0.8, d)
}
