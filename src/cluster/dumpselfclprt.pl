### perl ###
##!/usr/bin/perl -w

use dbinfo;
use Pg;
use Ace;

$pgdb = get_orparadb;

$in_table = "gpcr_cluster";
$query = "select rep, prt from $in_table order by rep, prt";

$ace_db = get_local_ace("human");

$result = $pgdb->exec($query);
if ($result->resultStatus != PGRES_TUPLES_OK) {
	die "$query failed on pg database\n";
}
@row = $result->fetchrow;
while (@row) {
	$rep = $row[0];
	open CLU, ">$rep.cl" or die "Cannot open $rep for write\n";

	print "\ncluster: $rep\n";
	while (@row && $rep eq $row[0]) {
		$prt = $row[1];
		$prt = $ace_db->fetch("Protein", $prt);
		$pep = $prt->asPeptide;
		$title = $prt->at("Title")->at;
		print "\t$prt\t$title\n";
		print CLU $pep;
		@row=$result->fetchrow;
	}
	close CLU;
}
