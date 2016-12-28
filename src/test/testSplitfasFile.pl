#!/usr/bin/perl -w

use Bioseq;

my ($infile, $binsize);
my $i=0;
while ($ARGV[$i]) {
    if ($ARGV[$i] eq "-s") { $binsize=$ARGV[++$i]; }
    else {
        $infile=$ARGV[$i];
    }
    ++$i;
}

splitFastaSeqfile($infile, $binsize);
