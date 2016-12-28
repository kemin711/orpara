### perl ###
##!/usr/bin/perl -w

$count=0;
$nw=0;
while (<>) {
   ++$nw;
	@row=split /\t/;
    if (! exists $ids{$row[0]}) {
      print ">", $row[0], "\n", $row[9], "\n";
      $ids{$row[0]}=1;
      ++$count;
    }
}

print STDERR "$nw rows $count sequences\n";
