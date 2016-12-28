package SeqIterator;

use Bioseq;
#use DBI;
use Exporter;
use strict;

our @ISA= ('Exporter');
our @EXPORT = qw(&validateFastaFile &checkDuplicatedSeqid);
my $DEBUG=1;

=head1 NAME

SeqIterator - Specialied to iterate through the fasta formated sequence files

=head1 SYNOPSIS

 my $reader=SeqIterator->new($nrfile, 1)
 if you want to collect organism information from the header.
    otherwise,
 my $reader=SeqIterator->new($nrfile, 0)

    $seq=$reader->next;
    $id=$reader->getPrimaryId;

=head1 Attributes

Normally the following attributes depends on the information
present in the fasta formated sequence header.  The sequence
could be

 >id
 HEQWWAACCEEAGMNH

or 

 >id title
 HEQWWAACCEEAGMNH
 LINEACT

 Data Structure

 headers is an array of strings storing potentially
 multiple headers. headers[0] is the primary header.

=head2 new

my $seqwalker = seqIterator->new($fasta_file, $orginfo);
This method will also load the species names from the 
taxonomy database which is stored locally as a mysql database.

$orginfo is a flag to indicate whether to load species
names from the database or not.
This is poor design, should getrid of it in the future.
I am getting rid of it.

Right now the taxonomy database is located in taxonomy@frylock
There is an identical copy on kemin@shake

=cut

sub new {
	my $invocant = shift;
    my $infile = shift;
    #my $orginfo=shift; # flag to get organism or not

	my $class = ref($invocant) || $invocant;
	my $fh;
	open($fh, "<", $infile) or die $!;
	my $self = { fh => $fh };
	my $buff = <$fh>; # first line
	if (!$buff) {
      die "Empty input file: $infile! inside seqIterator->new()\n";
    }
    elsif ($buff !~ /^>[A-Za-z0-9].*?\b/) { 
		die "not a fasta file: $buff inside seqIterator->new\n";
	}
# catch bad subheader division marks
    #my $badchar=ord(1); # ^A in vim program
    #my $tmp=$badchar . 'gi\|';
    #if ($buff =~ /$tmp/) {
    #    $self->{subhead_divider} = $badchar;
    #}
    #else {
    #    $self->{subhead_divider} = " >";
    #}
	chomp $buff;
	$self->{"nextHeader"} = $buff;
    $self->{infile} = $infile;
	bless($self, $class);
    #if ($orginfo) {
    #    $self->getOrgnames; # load species names from mysql DB: taxonomy
    #}
	return $self;
}

=head2 next

walks to the next sequence

return the next sequence in one string, not line terminated
This method does not check for illegal amino acids symbols.

Other informaton is accessed through other methods.

This function does not validate the sequence symbols
because it could be either amino acids or nucleic acids,
and the sequence text could use symbols other than the
known ones. The user has to do it by themselves.

=cut

sub next {
	my $self = shift;
	if (!($self->{"nextHeader"})) {
		return undef;
	}
    $self->processHeader($self->{nextHeader});

    #$self->{header} =~ s/\t+/ /g; # convert tab to space
# This is to fix a bug in certain fasta files. 
# there should be no tab in fasta files.
	my $fh = $self->{"fh"};

	my $seq = <$fh>;
	chomp $seq;
	my $ln = <$fh>;
	while ($ln && $ln !~ /^>/) {
		chomp $ln;
        if ($ln !~ /^\s*$/) { # make sure not empty line
            $seq .= $ln; 
        }
        $ln = <$fh>;
	}
	if (!$ln) { 
		#print STDERR "no more header\n";
		undef $self->{nextHeader};
	}
	else {
		chomp $ln;
		$self->{nextHeader} = $ln;
	}
	#$self->{header} =~ s/ >(\d)/ $1/;
	#$self->{header} =~ s/ > / /;
	#$self->digestHeader;
    $seq =~ s/ //g;
    $seq =~ s/[^A-Za-z]//g;
    $self->{seq}=$seq;
	return $seq;
}

sub U2T {
    my $self=shift;
    $self->{seq} =~ tr/U/T/;
}

=head2 header

return the first header line as text if there are multiple
headers such as NCBI nr file.

header has the > removed now 

=cut

sub header {
	my $self = shift;
    my $newh=shift;
    if ($newh) {
        $self->{headers}->[0]=$newh;
    }
    if (exists $self->{headers}->[0]) {
        return $self->{headers}->[0];
    }
    else { return ''; }
}

=head2 writeFas

write this object into a file handle
    The output have line terminator.

=cut
sub writeFas {
    my $self=shift;
    my $fh=shift;
    if (!$fh) { $fh=\*STDOUT; }
    print $fh $self->header, "\n";
    printSeq($self->{seq}, $fh, 70); 
    #print $fh "\n";
}

=head2 title

Get the first subheader without the > sign and without the id part
Some fasta files may not have titles!

Each fasta header may have more than one headers, especially
with the nr distribution where the header has been abused.

>gi|123|gb|GP09117.2| title1 junk words >gi|456|gb|AB0123| title2 words

This method will return the title1 words only if this exists

if the fasta file is the format

>id123

wihtout any title this method will return undef

gi|4966374|gb|AAD34705.1|AC006341_33
F3O9.28 [Arabidopsis thaliana]

got empty title!
Use of uninitialized value in print at /home/kzhou/proj64/seqana/seqsep/splitnr line 192, <$fh> line 13175735.


=cut
sub title {
    my $self=shift;
    if (exists $self->{headers}->[0] && 
        $self->{headers}->[0] =~ /^.+? (.+)$/) 
    {
        return $1;
        #return $self->{headers}->[0];
    }
    else {
        #$self->showHeaders;
        #print "\ngot empty title!\n";
    }
    return undef;
}

sub showHeaders {
    my $self=shift;
    print STDERR join("\n", @{$self->{headers}}), "\n";
}

=head2 subHeaders

return a reference to the subheader array
one header may have many >seqid1 header1 >seqid2 header2 .....
The '>' char is removed.

This method is specialized for the nr format
where multiple >seqid >seqidx are present

=cut

sub subHeaders {
	my $self=shift;
	#return $self->{subheaders};
    my @subh = @{$self->{headers}};
    if (@subh > 1) {
        shift @subh;
        return \@subh;
    }
	return undef;
}

=head2 getId

return the primary ID from the header of fasta formated files

The primary id is the first id encountered.

=cut

sub getId {
	my $self = shift;
	my $head = $self->header;
	if ($head =~ /gi\|(.+?)[| ]/ || $head =~ /IPI:(.+?)[| ]/ 
         || $head =~ /^>(.+?) / || $head =~ /^>(.+?)$/ 
         || $head =~ /^([^\s]+)$/  || $head =~ /^([^\s]+) /
            ) {
		return $1;
	}
	else {
		die "Cannot extract seqid from header: $head! inside getId()\n";
	}
}

=head2 getPrimaryId

returns the first id, in the following case the gi
>gi|28374049|pdb|1NLV|A

this will return the header.  In the case of NR, there are multiple
headers concatenated together.  This attribute store the first one.

=cut

sub getPrimaryId {
	my $self = shift;
	return $self->getId;
}

=head2 getAllIds

return an array of all the ids including the primary id,
the primary id is the first element in the array.

=cut

sub getAllIds {
	my $self=shift;
	my @ids=();
	foreach my $hd (@{$self->{headers}}) {
        #print STDERR "working on header: $hd ...\n";
        if ($hd =~ /^(.+?)\s/ || $hd =~ /^([^\s])$/) {
            my @iddbpair=digestIds($1);
            foreach my $id (@iddbpair) {
                push @ids, $id->[0];
            }
        }
        else {
            die "improper header $hd\n";
        }
	}
	return @ids;
}

=head2 processHeader($raw_header)

This is a helper function, and should not be used by the
external world.

This method does not do automatic processing of the headers.
It should be work on demand.

all subheadings got ">" removed.

{headers} => [ subheaders without > ]

=cut
sub processHeader {
    my $self=shift;
    my $rawheader=shift;
    $rawheader =~ s/\t+/ /g; # convert tab to space
    my $tmp=chr(1) . 'gi\|'; # probe for some bug in NCBI
    my @subh;
    if ($rawheader =~ /$tmp/) {
        $tmp=chr(1);
        @subh = split /$tmp/, $rawheader;
        #print STDERR "headers right after split:\n", join("\n", @subh), "\n\n";
    }
    else { @subh = split / >/, $rawheader; }
    $subh[0] =~ s/>//;
    $self->{headers}=\@subh;
}

=head2 getSpecies

    extract the first possible species inside [ ]

    Most of the time you only need one species name.
If there are multiples, then they are subspecies or strains.

=cut
sub getSpecies {
    my $self=shift;
    foreach my $header (@{$self->{headers}}) {
        #print STDERR "extracting species from header:\n$header ...\n";
        my $bi=rindex($header, '[');
        my $sp;
        if ($bi != -1) {
            #print STDERR "found [ at $bi\n";
            $sp=substr($header, $bi+1);
            $sp =~ s/\]$//;
        }
        elsif ($header =~ / - (.+?)$/) { $sp=$1; }
        else {
            #warn "header has no species: $header\n";
        }
        if ($sp && $sp =~ / - (.+?)$/) { $sp=$1; }
        if ($sp) { 
            $sp =~ s/ \(.+\)$//; 
            $sp =~ s/\s+$//;
        }

        if ($sp) {
            if (isSpecies($sp)) {
                return $sp;
            }
            else {
                #print STDERR "$sp is not species\n";
            }
        }
        else {
            #print STDERR "$header has no species\n";
        }
    }
    return undef;
}

=head2 getAllSpecies

return an array of all the species names.

the fasta file use [specie name] format
if [] is not present it will return an empty array

This methods is only usefull for nr database where
some of the entries have organism info appended
at the end [species name]

=cut

sub getAllSpecies {
	my $self=shift;
    my %species;
    my $numsp=0;
    foreach my $header (@{$self->{headers}}) {
        if ($header =~ /\[(.+)\]$/) {
            my $sp=$1;
            if (isSpecies($sp)) {
                ++$species{$sp};
                ++$numsp;
            }
        }
    }
    if ($numsp > 0) {
        my @species;
        foreach my $sp (keys %species) {
            push @species, $sp;
        }
        $self->{species}=\@species;
        return @species;
    }
	return ();
}

=head2 getAllIds_dbsrc 

return an array of [id, db] pairs
include the primary id as the first element.

=cut

sub getAllIds_dbsrc {
	my $self=shift;
    if (exists $self->{ids}) {
      return @{$self->{ids}};
    }
    else {
      return ();
    }
}


=head2 digestHeader

this method will be called automaticaly in the next method

This method will set up all the essential information
about headers.

Set up Attributes of
    ids
    species
    subheaders

=cut
sub digestHeader {
	my $self=shift;
# header always exits
    my $tmp=chr(1) . 'gi\|';
    my @subh;
    my %species=();
    my ($h, $ids, $des, $org);
    if ($self->{header} =~ /$tmp/) {
        $tmp=chr(1);
        @subh = split /$tmp/, $self->{header};
    }
    else {
        @subh = split / >/, $self->{header};
    }
    $h = shift @subh; # has removed it!
    ($ids, $des)=split(/\s+/, $h, 2);
	if ($des) {
		$self->{primaryheader}=$des;
		if ($org=extractSpecies($des)) {
			++$species{$org};
		}
	}
    my @iddb = digestIds($ids);
    if (@subh > 1) { # if there is subheaders
        $self->{subheaders}=[@subh];
    }
# subsequent subheaders 
	foreach my $sh (@subh) {
		($ids, $des)=split(/\s+/, $sh, 2);
		push @iddb, digestIds($ids);
		if (!$des) { next; }
		$org=extractSpecies($des);
		if (!$org) { next; }
		$species{$org}++;
	}
	$self->{ids}=\@iddb;
	my @sp=keys %species;
	$self->{species}=\@sp;
}



# helper function testing str is species or not
# this string does not have [] on either side.
sub isSpecies {
	my $str=shift;
	if (!$str) { return 0; }
    if ($str =~ /[\[\]]/) { return 0; }
	if ($str =~ /^[A-Z][a-z]+ [a-z]+$/ || $str =~ /^[A-Z][a-z]+ [a-z]+ /) {
		return 1;
	}
	if ($str =~ /^[A-Z][a-z]+ [a-z]+( \w+){1,4}$/ 
			|| $str =~ /^[A-Z][a-z]+ [a-z]+ [a-z]{2}\.( .+){1,3}$/
			|| $str =~ /^[A-Z][a-z]+ [a-z]+ (var|str)\. .+$/
			|| $str =~ /^[A-Z][a-z]+ [a-z]+( [.-:A-Za-z\d]+){1,2}$/
			|| $str =~ /.+?virus.*?$/i
			|| $str =~ /.*?uncultured .*$/i
			|| $str =~ /.*?unidentified .*?bacteri.*$/i
		) 
	{
		return 1;
	}
	elsif ($str =~ /aa$/ || $str =~ /\d+-\d+:\d+-\d+$/
		|| $str =~ /\d-[A-Za-z]+$/
		|| $str =~ /\)$/
		|| $str =~ /(protein|peptide) .+?$/i
		|| length($str) > 120 
	) {
		return 0;
	}
	elsif ($str =~ /^.\d/ 
			|| $str =~ /(protein|peptide|kinase)/i
			|| $str =~ /chain/i
			) {
            return 0;
    }
    elsif ($str =~ /sp\./
        || $str =~ /^[A-Z][a-z]+ [a-z]+ [^ ]+$/
        || $str =~ /^[A-Z][a-z]+ [a-z]+ (\w+){1,2}$/
        || $str =~ /Escherichia coli/
    ) {
        return 1;
    }
    else {
        return 0;
	}
}

=head digestIds

    gi|23335287|ref|ZP_00120524.1|
    gi|13878750|sp|Q9CDN0.1|RS18_LACLA

    This method is not very easily implemented because 
    sp can have more than one ids!

    produce a [id, db] pair array

=cut
sub digestIds {
	my $str=shift;
	$str =~ s/\|\|/\|/g;  # replace double || with |
	$str =~ s/\|$//; # removal trailing |
	my @arr=split /\|/, $str;
    my $i=0;
    my @ids;
    my $db;
    while ($i < @arr) {
        if ($arr[$i] eq 'gi' || $arr[$i] eq 'gb' || $arr[$i] eq 'ref'
            || $arr[$i] eq 'sp' || $arr[$i] =~ /^[a-z]{2,4}$/)
        {
            $db=$arr[$i++];
        }
        push @ids, [$arr[$i], $db];
        ++$i;
    }
	return @ids;
}

=head2 extractSpecies

extract species name from the input string.  It first tries
to use regular expression to extract string contained inside
[ xxxx ]
Then it tries to figure out whether it is in the format of
species name or not. 

The string from [] may not be real species name stored in 
the taxonomy database! We are not checking it.

=cut

sub extractSpecies {
	my $str=shift;
	if (!$str) {
		warn "header is empty returning undef\n";
		return undef;
	}
# must appear at the end of the title
	if ($str =~ /\[([A-Z][a-z]+ [a-z]+)\]$/) {
		return $1;
	}
	if ($str =~ /\[([A-Z][a-z]+ [a-z]+( \w+){1,4})\]$/ 
			|| $str =~ /\[([A-Z][a-z]+ [a-z]+ [a-z]{2}\.( .+){1,3})\]$/
			|| $str =~ /\[([A-Z][a-z]+ [a-z]+ (var|str)\. .+)\]$/
			|| $str =~ /\[([A-Z][a-z]+ [a-z]+( [.-:A-Za-z\d]+){1,2})\]$/
			|| $str =~ /\[(.+?virus.*?)\]$/i
			|| $str =~ /\[(.*?uncultured .*)\]$/i
			|| $str =~ /\[(.*?unidentified .*?bacteri.*)\]$/i
		) 
	{
		#print STDERR $1, " ****\n";
		return $1;
	}
	elsif ($str =~ /aa\]$/ || $str =~ /\d+-\d+:\d+-\d+\]$/
		|| $str =~ /\d-[A-Za-z]+\]$/
		|| $str =~ /\)\]$/
		|| $str =~ /(protein|peptide) .+?\]$/i
		|| length($str) > 120 
	) {
		return undef;
	}
	elsif ($str =~ /\[(.+?)\]$/) {
		$str=$1;
		if ($str =~ /^.\d/ || $str =~ /\[/
			|| $str =~ /(protein|peptide|kinase)/i
			|| $str =~ /chain/i
			) {
# ignore
		}
		elsif ($str =~ /sp\./
			|| $str =~ /^[A-Z][a-z]+ [a-z]+ [^ ]+$/
			|| $str =~ /^[A-Z][a-z]+ [a-z]+ (\w+){1,2}$/
			|| $str =~ /Escherichia coli/
		) {
			return $str;
		}
		#elsif (exists $self->{spname}{$str}) {
		#	return $str;
		#}
		#else {
		#	print "title $str has no organism, code renew\n";
		#}
	}
	return undef;
}

=head2 getOrgnames

It reaches to a mysql database that has the taxonomy database schema.
Grab all organism names.  This is used for parsing the headers of the 
NCBI nr file because it contains the organism names.  This is a very
expensive operation.

This method should not be here, and I am removing it


sub getOrgnames {
	my $self=shift;
	my $dbh=DBI->connect("dbi:mysql:kemin:shake", "internal", "w00kie", {PrintError => 0, AutoCommit => 1} ) or die DBI::errstr;
	my $sql="select m.name from tax_nodes n join tax_names m on n.id=m.id where n.rank like '%species%'";
	my $sth=$dbh->prepare($sql);
	$sth->execute;
	my $count=0;
	my %spname=();
	while (my $rowref=$sth->fetchrow_arrayref) {
		$spname{$rowref->[0]}++;
		$count++;
	}
	print STDERR "$count species names loaded for look up\n";
	$self->{spname}=\%spname;
}

=cut

sub getBriefHeader {
	my $self=shift;
	my $hs = $self->subHeaders;
	my $numseq = scalar(@$hs);
	my @orgs=$self->getAllSpecies;
	my $brief=">" . $hs->[0] . "; Representing $numseq sequences";
	if (@orgs > 0) {
		$brief .= " from species: " . join(', ', @orgs);
	}
	return $brief;
}

# input reference of arrays
sub fetchSubset {
    my $self = shift;
    my $seqids = shift;
    my $outfh = shift;
    my %h = ();
    my $id;
    foreach $id (@$seqids) {
        $h{$id} = 1;
    }
    while (my $seq = $self->next) {
        $id = $self->getId;
        if (exists $h{$id}) {
            print $outfh $self->header, "\n";
            printSeq($seq, $outfh);
            delete $h{$id};
        }
        #print scalar( keys %h ), " seq to fetch ..\n";

        if (scalar(keys %h) < 1) {
            last;
        }
    }
    if (scalar(keys %h) > 0) {
        foreach my $k (keys %h) {
            print STDERR $k, "\n";
        }
        print STDERR " were not found in the input file: ", 
            $self->{infile}, "\n";
    }
}

sub DESTROY {
	my $self = shift;
	my $fh = $self->{"fh"};
	close $fh;
}

=head2 validateFastaFile

to make sure that the sequence does not have
inalid symbols [^A-Za-z] meaning not letters.

Right now only validating sequence, by removing
none letter characters, usually white char.

return the file name that is good. It could be
itself or a new file.

=cut
sub validateFastaFile {
    my $file=shift;
    open IN, "<$file" or die $!;
# remove foreign line terminator
    open TMP, ">tmp.fas" or die $!;
    my $foreignTerm=0;
    while (<IN>) {
        chomp;
        if (s/\015//g) { ++$foreignTerm; }
        s/\s+$//;
        print TMP $_, "\n";
    }
    if ($foreignTerm) {
        print STDERR "$foreignTerm foreign terminators removed\n";
        system("mv tmp.fas $file");
    }
    else { unlink 'tmp.fas'; }
###################################
##### check for non-aa characters #######

    my $goodfile=$file;
    $goodfile =~ s/\..+$//;
    $goodfile .= "_N.fas";
    open OU, ">$goodfile" or die $!;
    my $badcnt=0;
    my $logfile="invalidSeq.log";
    open LOG, ">$logfile" or die $!;

    my $it=SeqIterator->new($file);
    my $pepseq=$it->next;
    while ($pepseq) {
        $pepseq=cleanString($pepseq);
        $it->header(cleanString($it->header));
        if ($pepseq =~ s/([^A-Za-z])//g) {
            if ($DEBUG) {
                print LOG "non-amino acid symbol: '$1' found in sequence ",
                    $it->getId, "  $pepseq\n";
            }
            ++$badcnt;
        }
        print OU $it->header, "\n";
        printSeq($pepseq, *OU);
        $pepseq=$it->next;
    }
    if ($badcnt > 0) { 
        print STDERR "original file $file have $badcnt bad sequences, using a cleaned up version $goodfile\nInvalid sequences written into $logfile\n";
        return $goodfile;
        #return 0; 
    } # bad
    else { 
        unlink $goodfile;
        unlink $logfile;
        return $file; 
    }
}

sub cleanString {
    my $str=shift;
    $str =~ s/[[:cntrl:]]//g;
    $str =~ s/\cM//g;
    $str =~ s/\015//g;
    return $str;
}

=head2 checkDuplicatedSeqid

    @param $fasta_file, $output_fh
    The output_fh is optioal. If not given it will
    write to STDERR.

    @return the number of duplicated ids

=cut
sub checkDuplicatedSeqid {
    my $fasfile=shift;
    my $fh=shift;
    if (!$fh) { $fh=\*STDERR; }
    my @header=`grep '>' $fasfile`;
    chomp @header;
    my %ids;
    foreach my $h (@header) {
        if ($h =~ /^>(.+?) .+$/ || $h =~ /^>(.+)$/) {
            $ids{$1}++;
        }
        else {
            die "invalide fasta header $h\n";
        }
    }
    my $dupcnt=0;
    foreach my $i (keys %ids) {
        if ($ids{$i} > 1) {
            print $fh $i, " is duplicated ", $ids{$i}, " times\n";
            ++$dupcnt;
        }
    }
    return $dupcnt;
}

1;
