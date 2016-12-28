### perl ###
##!/usr/bin/perl -w

# mark the gi divisions using the following 
# Metazoa (Animals), Plant, Fungi, Eukaryota, Archaea, Bacteria

use Taxonomy;
use integer;
use strict;

my $taxdatadir="/house/groupdirs/gat/projects/kzhou/taxonomy";
loadTaxonomy($taxdatadir);
# markdiv();
# now moved into Taxonomy package, you don't have to 
# do this manually.
print "Done!\n";

__END__

please refere to proj/automation/loadtaxonomy for
loading programs, that loads the results of download phase
which is accomplished through manual work.  My automation
is having some file format problems.
