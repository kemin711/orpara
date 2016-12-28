#!/usr/bin/perl -w

use SeqIterator;
use Bioseq;

#test_validateFastaFile();
test_iterator("/house/blast_db/nrJuly2010/nrhead50.fas");
print "Done!\n";

###################################

sub test_iterator {
    my $fasfile=shift;
    my $it=SeqIterator->new($fasfile);

    my $seq=$it->next;
    for (my $i=0; $i<10; $i++) {
        print "primary Id: ", $it->getId, "\n"; 
            #"all ids: ", join(' | ', $it->getAllIds), "\n";
        print "species: ", $it->getSpecies, "\n";
        $seq=$it->next;
    }
}

sub test_validateFastaFile {
    print STDERR "testing validateFastaFile ...\n";
    my $file="/home/kzhou/work/pepmap/Sporotrichum_thermophile/frame/St_pep6frame.fas";
    my $validfas=validateFastaFile($file);
    system("cp $validfas .");
}
