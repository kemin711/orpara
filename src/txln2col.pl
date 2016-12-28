#!/usr/bin/perl -w



my $i=0;
while ($ARGV[0]) {
   my $tmp=shift;
   $tmp =~ s/,\s*//;
   push @arr, lc($tmp);
}
print "\n\n", join("\n", @arr), "\n";

__END__

=head1 DESCRIPTION

this simple utility converts text with comma separated fields
into a column data:

GENOME_DNA_ID, SOURCE_ID, CHROMOSOME_ID, GI_NO, ACCESSION_NO,

to

GENOME_DNA_ID
SOURCE_ID
CHROMOSOME_ID
GI_NO
ACCESSION_NO

I designed this application to deal with Tode could not copy column
names.
