#!/usr/bin/perl -w

open IN1, "<testcpp.ltrace.out" or die $!;
open IN2, "<yy" or die $!; # crashed testcpp with ltrace -d testseq

my $i=1;
my $j=1;

my ($ln1, $ln2);

$ln1=<IN1>;
$ln2=<IN2>;
while ($ln1 && $ln2 && $i < 508 && $j < 508) {
    chomp $ln1;
    chomp $ln2;
    $ln1 =~ s/^DEBUG: //;
    $ln2 =~ s/^DEBUG: //;
    if ($ln1 eq $ln2) {
        $ln1=<IN1>;
        $ln2=<IN2>;
        ++$i; ++$j;
    }
    else {
        my $prefix = commonPrefix($ln1, $ln2);
        if (length($prefix)/length($ln1) > 0.5 
                && length($prefix)/length($ln2) > 0.5) {
            print "$i |$ln1|\n", "$j |$ln2|\n", "     Different prefix\n";
            $ln1=<IN1>;
            $ln2=<IN2>;
            ++$i; ++$j;
        }
        else {
            print "Quite different:\n$ln1\n$ln2\n";
            exit(0);
        }
    }

}

sub commonPrefix {
    my $l1=shift;
    my $l2=shift;
    my @arr1 = split //, $l1;
    my @arr2 = split //, $l2;

    my $i=0;
    while ($i < @arr1 && $i < @arr2) {
        if ($arr1[$i] eq $arr2[$i]) {
            ++$i;
        }
        else {
            last;
        }
    }
    return substr($l1, 0, $i);
}

