### perl ###
##!/usr/bin/perl -w

# given a number of files this program
# will map the id to numbers
# in the end it will generate a table stringid => integerid

use strict;

my @files;
my $counter=1;
my $outfile='idindex.tab';

my $i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq '-s') { $counter=$ARGV[++$i]; }
	elsif ($ARGV[$i] eq '-o') { $outfile=$ARGV[++$i]; }
	else {
		push @files, $ARGV[$i];
	}
	++$i;
}

open OU, ">$outfile" or die "Failed to open output file $outfile $!\n";
foreach my $f (@files) {
    open IN, "<$f" or die $!;
    print STDERR "Working on $f ...\n";
    while (<IN>) {
        if (/^>([A-Za-z0-9_]+)\b/) {
            print OU "$1\t", $counter++, "\n";
        }
    }
}
print STDERR "Index written to $outfile lastid $counter\n";

