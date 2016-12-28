### perl ###
###!/usr/bin/perl -w

use SeqIterator;
use Statistics::Welford;
use strict;

my $i=0;
my @files;
my $lendistFile='seqlength.tab';
while ($ARGV[$i]) {
   if ($ARGV[$i] eq '-o') { $lendistFile = $ARGV[++$i]; }
   else {
      push @files, $ARGV[$i];
   }
   ++$i;
}
if (@files<1) {
    usage();
}

# global will accumulate stat from all given files
my $stat = Statistics::Welford->new;
my %lenstat;
#open LEN, ">length.tab" or die $!;

foreach my $f (@files) {
    processOneFile($f);
}
#close LEN;

## for plotting
close OU;
open OU, ">$lendistFile" or die $!;
foreach my $k (sort {$a <=> $b } keys %lenstat) {
    print OU $k, "\t", $lenstat{$k}, "\n";
}
print STDERR "Length distribution written to $lendistFile\n";
print "min\tmean\tmax\tstd\tcount\n",
    $stat->min, "\t", $stat->mean, "\t", $stat->max, "\t", 
    $stat->standard_deviation, "\t", $stat->n, "\n";

############# Subroutines ###############

sub usage {
    print STDERR "Usage faseqlenstat -o lengthDistributionFile <fasfile1> <fasfile2> ...\n",
         "If the lengthDistributionFile name is not given then it will use the default\n";
    exit(1);
}

sub processOneFile {
    my $fasfile=shift;
    my $reader=SeqIterator->new($fasfile);
    my $len;

    while (my $seq=$reader->next) {
        $len=length($seq);
        #print LEN $reader->getId, "\t", $len, "\n";
        ++$lenstat{$len};
        $stat->add($len);
    }
}


__END__

=head1 NAME

faseqlenstat - comput simple statistics about sequence length

=head1 SYNOPSIS

faseqlenstat seq1.fas seq2.fas ...

=head1 DESCRIPTION

Given a few fasta files, this program will produce statistics for the 
length of the sequences.
