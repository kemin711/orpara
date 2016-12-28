#!/usr/local/bin/perl -w

use strict;
use Basicstat;

my @values = (1.3, 7.5, 8.4, 2.7, 4.4, 3.2);
my @stat = avgstd_pop(\@values);
my @returnType = getStatValueType;
print join(' | ', @returnType), "\n";
print join(' | ', @stat), "\n";

@values = (2.7, 4.4);
@stat = avgstd_pop(\@values);
@returnType = getStatValueType;
print join(' | ', @returnType), "\n";
print join(' | ', @stat), "\n";
