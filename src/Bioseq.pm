package Bioseq;

# trying to convert the package into an
# object orientated package, the procedure orientated
# methods will also stay, I will also add some
# methods that are more appropriate for simple methods

use Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(&getSeqFromFile &nameOutf &revcomp &translate &convertFas2String &addFasHeader &printSeq &splitFastaSeqfile &complement &getFasPrimaryId &getFasAllId &getFasAllId_dbsrc readFastaSeq &subseq &subsequence &fileStem &breakfasta &setCodonTable);

# indexed by the codon table id [codon hash codon=>AA, start codons]
my @codontables;
my $codontable=1; # default universal


# read codon table
BEGIN {
    #our $codontableFile = $ENV{HOME} . "/home/kzhou/etc/codontable.txt";
    #our $codontableFile = $ENV{HOME} . "/etc/codontable.txt";
   our $codontableFile = "/remote/RSU/sw-cache/metag/lib/perl5/site_perl/codontable.txt";
   # load codon tabl from file
   open CIN, "<$codontableFile" or die "Cannot open codon table file $codontableFile\n";
   my $ln=<CIN>;
   while ($ln) {
      chomp $ln;
      if ($ln =~ /^\s*$/) { $ln=<CIN>; }
      elsif ($ln =~ m|^//|) { $ln=<CIN>; }
      elsif ($ln =~ m|^/\*|) { 
         if ($ln !~ m|\*/$|) {
            $ln=<CIN>; 
            while ($ln && $ln !~ m|\*/|) {
               chomp $ln;
               $ln = <CIN>; 
            }
         }
         $ln=<CIN>;
      }
      else {
         my $ctid; # codon table id
         my @bases=('T', 'C', 'A', 'G');
         if ($ln =~ /^(\d+)\. /) { $ctid=$1; }
         else {
            die "must be codon table header line\n";
         }
         #print "working on $ctid ...\n";
         $ln=<CIN>;
         while ($ln !~ /^  AAs  = /) { $ln=<CIN>; }
         chomp $ln; $ln=substr($ln, 9);
         my $startln=<CIN>;
         chomp $startln; $startln=substr($startln, 9);
         my @pep=split //, $ln;
         my @sss=split //, $startln;
         my $x=0;
         my (%codon, %starts);
         for (my $i=0; $i<4; ++$i) {
            for (my $j=0; $j<4; ++$j) {
               for (my $k=0; $k<4; ++$k) {
                  my $codtxt=$bases[$i] . $bases[$j] . $bases[$k];
                  $codon{$bases[$i] . $bases[$j] . $bases[$k]}=$pep[$x];
                  if ($sss[$x] eq 'M') {
                     $starts{$codtxt}=1;
                  }
                  ++$x;
               }
            }
         }
         $codontables[$ctid]=[\%codon, \%starts];
         $ln=<CIN>;
         while ($ln && $ln =~ /^Base[123]  = /) {
            $ln=<CIN>;
         }
      }
   }
   #print "codon table in initialization block\n";
   #printValidCodonTables();
}

sub printValidCodonTables {
   print scalar(@codontables), " Valid codons tables\n";
   for (my $i=0; $i < @codontables; ++$i) {
      if ($codontables[$i]) {
         print "$i  ";
      }
   }
   print "\n";
}

=head2 getCodonTable

obtain a particular codon table by id

=cut
sub getCodonTable {
   my $cid=shift;
   if (!$codontables[$cid]) {
      print $codontables[$cid], "\n";
      die "codon table $cid does not exist\n";
   }
   return $codontables[$cid]->[0];
}

=head2 setCodonTable

   set the global codon table id 

=cut
sub setCodonTable {
   #print "codon table inside setCodonTable()\n";
   #printValidCodonTables();
   print "setting codon table to ", $_[0], "\n";
   $codontable=$_[0];
}

=head1 NAME

bioseq - A collection of convenient function dealing with biosequences.

=head1 SYNOPSIS

splitFastaSeqfile($infile, $numseqPerFile)

=head1 DESCRIPTION

Here each method will be described clearly.

=head2 readFastaSeq($fh, $terminator, $nextline) 

Read one fasta sequence from the file stream
using terminator to tell the end of the input.

It can start with a header 
>seq1
FDAGK....

If there are multiple sequences, then it will
pass out the next line in the third argument.

return ($title,$seqstr)

=cut

sub readFastaSeq { 
    my $fh=$_[0];
    my $term=$_[1];

    my $seq = "";
    my $title = <$fh>;
    chomp $title;
    if ($title !~ /^>/) {
        # the title is missing
        $seq .= $title;
        $title="";
    }

    my $l = <$fh>;
    while ($l && $l !~ /^($term|>)/o) {
        chomp $l;
        #$l =~ s/\W+//g;  #pick sequence only avoiding ^M thing
        $l =~ s/[^A-Za-z*]+//g;  #pick sequence only avoiding ^M thing
        $seq .= $l;
        $l = <$fh>;
    }
    chomp $l;
    $_[2]=$l;
    return ($title, $seq);
}

#sub new {
#    my $invocant = shift;
#    my $class = ref($invocant) || $invocant;
#    my $self = { @_ };
#    bless($self, $class);
#    return $self;
#}
#sub seqIterator {
#    # file name as argument
#    my $self = $_[0];
#    my $fh;
#    open($fh, $_[1]);
#    $self->{"fh"} = $fh;
#}

=head2 getSeqFromFile($fileName)

return a sequence from a file default fasta file only
the sequence is in a string without space or \n
calling getSeqFromFile([file name], [string of first line])
only ONE sequence per File!!!!!

=cut

sub getSeqFromFile { 
    my $infile = $_[0];
    open(INF, $infile) or die "Can't opn $infile: $!\n";
    my $l = <INF>;
    chomp $l;
    if ($l =~ /^>/) {
        $_[1] = $l;  # pass out the title of the sequence, first line
    }
    else {
        die "sequence not in fasta format!\n$l\n";
    }
    #assumes that the first line is >seqname etc
    my $seq = "";
    $l = <INF>;
    while ($l && $l !~ /^>/) {  #only reads the first sequence
        #chomp($l);
        #if (/-/) {
        #    die "- detected in sequence\n";
        #}
        $l =~ s/[^A-Za-z*]+//g;  #pick sequence only avoiding ^M thing
        $seq .= $l;
        $l = <INF>;
    }
    return $seq;
}

=head2 nameOutf($infile, $extension)

name the output file, find the last period
xxx.yy.ext => xxx.yy.new_ext

This is a helper function to manipulate file names.
It is commonly used function in bioinformatics.

=cut

sub nameOutf {
#(infile, extension), returns outfile with .extension

    my $infile = $_[0];
    my $new_ext = $_[1];

    if (!($infile =~ s/\.\w+?$/\.$new_ext/)) {
        $infile .= ".$new_ext";
    }
    return $infile;
}

sub complement {   
#return a new sequence, the original sequence is intact
    #use less memory
    #used for char based algorithms
   #my %comp= ('A', 'T', 'C', 'G', 'G', 'C', 'T', 'A', 'R', 'Y', 
    #    'Y', 'R', 'K', 'M', 'M', 'K', 'S', 'S', 'W', 'W', 'B', 'V',
    #   'D', 'H', 'H', 'D', 'V', 'B', 'N', 'N',
    #      'a', 't', 'c', 'g', 'g', 'c', 't', 'a', 'r', 'y', 'y', 'r',
    #      'k', 'm', 'm', 'k', 's', 's', 'w', 'w', 'b', 'v',
    #      'd', 'h', 'h', 'd', 'v', 'b', 'n', 'n');

    #decemal representation, work with both lower and upper case
   my %compnum= (65, 84, 67, 71, 71, 67, 84, 65, 82, 89, 89, 82, 
        75, 77, 77, 75, 83, 83, 87, 87, 66, 86, 68, 72, 72, 68, 86, 
        66, 78, 78,
       97, 116, 99, 103, 103, 99, 116, 97, 114, 121, 121, 114, 107,
       109, 109, 107, 115, 115, 119, 119, 98, 118, 100, 104, 104, 
       100, 118, 98, 110, 110);
   
   my $cseq = $_[0];
    if (!$cseq) { #empty string
        return "";
    }
    #my $len = length($rcseq);    
    for (my $i=0; $i<length($cseq); $i++) {
      vec($cseq, $i, 8) = $compnum{vec($cseq, $i, 8)};
    }
   return $cseq; 
}

=head2 revcomp

the sequence must be just a string

usage revcomp($original_String)

=cut
sub revcomp {   
#return a new sequence, the original sequence is intact
    #use less memory
    #used for char based algorithms
   #my %comp= ('A', 'T', 'C', 'G', 'G', 'C', 'T', 'A', 'R', 'Y', 
    #    'Y', 'R', 'K', 'M', 'M', 'K', 'S', 'S', 'W', 'W', 'B', 'V',
    #   'D', 'H', 'H', 'D', 'V', 'B', 'N', 'N',
    #      'a', 't', 'c', 'g', 'g', 'c', 't', 'a', 'r', 'y', 'y', 'r',
    #      'k', 'm', 'm', 'k', 's', 's', 'w', 'w', 'b', 'v',
    #      'd', 'h', 'h', 'd', 'v', 'b', 'n', 'n');

    #decemal representation, work with both lower and upper case
   my %compnum= (65, 84, 67, 71, 71, 67, 84, 65, 82, 89, 89, 82, 
        75, 77, 77, 75, 83, 83, 87, 87, 66, 86, 68, 72, 72, 68, 86, 
        66, 78, 78,
       97, 116, 99, 103, 103, 99, 116, 97, 114, 121, 121, 114, 107,
       109, 109, 107, 115, 115, 119, 119, 98, 118, 100, 104, 104, 
       100, 118, 98, 110, 110);
   
   my $rcseq = $_[0];
    if (!$rcseq) { #empty string
        return "";
    }
    #my $len = length($rcseq);    
   my $j=length($rcseq)-1;
   my $i=0;
    my $bc;
   while ($i < $j) { 
      $bc = $compnum{vec($rcseq, $i, 8)};
        vec($rcseq, $i, 8) = $compnum{vec($rcseq, $j, 8)};
        vec($rcseq, $j, 8) = $bc;
        $i++;
        $j--;
   }
    if ($i == $j) {  #when there are odd number of bases
        vec($rcseq, $i, 8) = $compnum{vec($rcseq, $i, 8)};
    }
   return $rcseq; 
}
######## not using array, may break in the future
## another straitforward implementation is as follows
#
#   %comp= ('A', 'T', 'C', 'G', 'G', 'C', 'T', 'A', 'R', 'Y', 'Y', 'R',
#            'K', 'M', 'M', 'K', 'S', 'S', 'W', 'W', 'B', 'V',
#             'D', 'H', 'H', 'D', 'V', 'B', 'N', 'N');
#    #$i=length($seq)-1;
#    @seqarr=split //, $seq;
#    $i=$#seqarr;
#    $j=0;
#    while ($i > -1) {
#       $rcseq[$j++]=$comp{$seqarr[$i]};
#       $i--;
#    }
#    $seq = join '', @rcseq;
#}

=head2 printSeq($seq, *fileHandle, $width)

    Input: $seq as a single string, 
           $width is the number of residues per line.
    Output: fileHandle

    70 characters per line in output
    otherwise can be spcified in the third argument

=cut
sub printSeq {
    my $seq = $_[0];
    my $fh=$_[1];
    if (!$fh) {
        $fh=\*STDOUT;
    }

    my $i=0;
    my $len = $_[2] ? $_[2] : 70;
    my $sub = substr($seq, $i, $len);
    my $seqLen = length($seq);
    print $fh "$sub\n";
    $i += $len;
    while ($i < $seqLen) {
        $sub = substr($seq, $i, $len);
        print $fh "$sub\n";
        $i += $len;
    }
}

sub addFasHeader {
#usage addFasHeader($seq, "additional text")
#add more text info to the header of fasta file
    my $fasseq = $_[0];
    my $i = index($fasseq, "\n");
    substr($fasseq, $i, 0) = " " . $_[1];
    return $fasseq;
}

# returns a string the is the sequence body
sub convertFas2String {
    #remove non-base character, and the header
    #returns the Pure DNA or Peptide as a long string
    #the header is stored in the second argument if provided
    # otherwise it will be discarded
    my $str = $_[0];
    $str =~ s/(^>.+?)\n//;
    if (@_ > 1) { $_[1] = $1; }
    $str =~ s/\W+//g;   # remove non word char
    return $str;
}

=head2 translate($dna, $start, $end)

translate DNA sequence into protein
usage $pep = translate($seq, $start, $end)

the input sequence must be string,
$atart and $end are optional and use 1-based index.

The implementation is a lot faster than the more generic version using
translation table. This method will use the codon table method if the codon
table is set to be other than 1 (the default codon, or universal codon table).

=cut
sub translate {
    my $s = $_[0];
   if (!$s) {
        die "empty input DNA sequence!\n";
   }
   my $start = $_[1];
    my $end = $_[2];
   # this method is only defined for the universal codon table.
   if ($codontable != 1) {
      return translate_codontable($s, $start, $end);
   }
    #      codon table state machine
   #
   #                                 Codon Table
   #
   #            First U       C       A            G            Last
   #
   #            U     Phe     Ser     Tyr          Cys          U
   #
   #                  Phe     Ser     Tyr          Cys          C
   #
   #                  Leu     Ser     Stop (Ochre) Stop (Umber) A
   #
   #                  Leu     Ser     Stop (Amber) Trp          G
   #
   #            C     Leu     Pro     His          Arg          U
   #
   #                  Leu     Pro     His          Arg          C
   #
   #                  Leu     Pro     Gln          Arg          A
   #
   #                  Leu     Pro     Gln          Arg          G
   #
   #            A     Ile     Thr     Asn          Ser          U
   #
   #                  Ile     Thr     Asn          Ser          C
   #
   #                  Ile     Thr     Lys          Arg          A
   #
   #                  Met     Thr     Lys          Arg          G
   #
   #            G     Val     Ala     Asp          Gly          U
   #
   #                  Val     Ala     Asp          Gly          C
   #
   #                  Val     Ala     Glu          Gly          A
   #
   #                  Val     Ala     Glu          Gly          G
   #
   #           ------------------------------------------------------
   #                                Chris Seidel
   #
   #                                    Home
    #  
    #      Start State | next State
    #      --------------------------------------------------
    #      ----------Input Character
    #        0   U     C      A    G
    #  ============================
    #          0   1     2      3    4
    #          1   5     6++S   7    8
    #          2   9++L  10++P  11   12++R
    #          3   13    14++T  15   16
    #          4   17++V 18++A  19   20++G
    #          5   21F   21F    22L  22L
    #          7   23Y   23Y    24*  24*
    #          8   25C   25C    24*  26W
    #          11  27H   27H    28Q  28Q
    #          13  29I   29I    29I  30M
    #          15  31N   31N    32K  32K
    #          16  33S   33S    34R  34R
    #          19  35D   35D    36E  36E
    #  
    #        i   U     C      A    G     N     Y(C+T)   R(A+G) 
    #  =========================================================
    #          0   1     2      3    4    ++X    ++X      ++X
    #          1   5     +S     6    7    +X      +X      +X
    #          2   +L    +P     8    +R   +X      +X      +X
    #          3   9     +T     10   11   +X      +X      +X
    #          4   +V    +A     12   +G   +X      +X      +X
    #          5   F     F      L    L      X     F       L
    #          6   Y     Y      *    *      X     Y       *
    #          7   C     C      *    W      X     C       X
    #          8   H     H      Q    Q      X     H       Q
    #          9   I     I      I    M      X     I       X
    #          10  N     N      K    K      X     N       K
    #          11  S     S      R    R      X     S       R
    #          12  D     D      E    E      X     D       E

    # + states indicate skip the next nucleotide

    my @stateArr = (
         [1,     2,     3,    4,    "++X", "++X", "++X"],
         [5,     "+S",  6,    7,    "+X",  "+X",  "+X"],
         ["+L",  "+P",  8,    "+R", "+X",  "+X",  "+X"],
         [9,     "+T",  10,   11,   "+X",  "+X",  "+X"],
         ["+V",  "+A",  12,   "+G", "+X",  "+X",  "+X"],
         ["F",   "F",   "L",  "L",  "X",   "F",   "L"],
         ["Y",   "Y",   "*",  "*",  "X",   "Y",   "*"],
         ["C",   "C",   "*",  "W",  "X",   "C",   "X"],
         ["H",   "H",   "Q",  "Q",  "X",   "H",   "Q"],
         ["I",   "I",   "I",  "M",  "X",   "I",   "X"],
         ["N",   "N",   "K",  "K",  "X",   "N",   "K"],
         ["S",   "S",   "R",  "R",  "X",   "S",   "R"],
         ["D",   "D",   "E",  "E",  "X",   "D",   "E"]   );

    my %ctoi = ( U=>0, u=>0, T=>0, t=>0,
                 C=>1, c=>1,
                             A=>2, a=>2,
                             G=>3, g=>3,
                             N=>4, n=>4, K=>4, k=>4, M=>4, m=>4, W=>4, w=>4, 
                             S=>4, s=>4, V=>4, v=>4, H=>4, h=>4, D=>4, d=>4,
                             B=>4, b=>4, O=>4, o=>4, Q=>4, q=>4, Z=>4, z=>4,
                             Y=>5, y=>5,
                             R=>6, r=>6    );

    my @seqArr = split //, $s;
    my $i = 0;
    my $input;
    if ($start) { $i = $start - 1; }  # 1 based to 0 based 
    if (!$end || $end > @seqArr) { $end = @seqArr; }
    my $state = 0;
    my $pep = "";
    while ($i < $end) {
        $input = $ctoi{$seqArr[$i]};  # input number
        if ($input !~ /[0-6]/) {
         #print STDERR "pos: $i base:",  $seqArr[$i], " Not in ctoi hash table\n";
            # print STDERR "sequence: $s\n";
            # print STDERR %ctoi, "\n";
            return 0;
            # die;
        }
        # $state = $stateArr[$state][$ctoi{$seqArr[$i]}];
        $state = $stateArr[$state][$input];   # for debug version
        if ($state !~ /[0-9]+/) {          # generating output
            if ($state =~ s/\+\+//) { $i += 2; }
            elsif ($state =~ s/\+//) { $i++; }
            $pep .= $state;
            $state = 0;
        }
        $i++;
    }
    return $pep;
}

sub translate_codontable {
   my $seq=shift;
   my $start=shift;
   my $end;
   if (!$start) { 
      $start=1;
      $end=length($seq);
   }
   else { 
      $end=shift;
      if (!$end) { $end=length($seq); }
   }

# use this method if codon table is not 1
   #if ($codontable == 1) {
   #   die "should have used the fast tranlate function!\n";
   #}
   my $ct=getCodonTable($codontable);
   my $i=$start-1;
   my $pep;
   while ($i+2 < $end) {
      if (! exists $ct->{substr($seq, $i, 3)}) {
         $pep .= 'X';
         #print "codon ", substr($seq, $i, 3), " not in codon table\n";
      }
      else {
         $pep .= $ct->{substr($seq, $i, 3)};
      }
      $i += 3;
   }
   #print "i $i  length(seq) ", length($seq), "\n";
   if (length($seq) != $i) {
      $pep .= "?";
   }
   return $pep;
}

=head2 subseq($seq, $ex) 

$ex is a reference to an array of exons location start with 1
or 0, this function will figure out this

for example: 
    1   50
    100 200
    500 570
 
=cut
sub subseq {
    my $seq=shift;
    my $ex=shift; # reference to an array of [start, end]
    if (length($seq) < 2) {
        die "length of sequence in Bioseq::subseq is zero!\n";
    }
    my $index_start=0;
    if ($ex->[0][0] == 1) {
        $index_start=1;
    }
    elsif ($ex->[0][0] > 1) {
        die "expecting 1-based index, first value must be 1\n";
    }
    my $buff='';
    foreach my $x (@$ex) {
        if ($x->[1] > length($seq)) {
            die "end of position, ", $x->[1], 
                " out of length of the source sequence whose lenghth is ",
                length($seq), " inside Bioseq::subseq\n";
        }
        $buff .= substr($seq, $x->[0]-$index_start, $x->[1]-$x->[0]+1);
    }
    return $buff;
}

=head2 subsequence

    input: (\$seq, \@exon)
    return \$subseq

    exon is an array of array [b1, e1], [b2, e2], ....
    The exons can also be given as a string in the format
        b1-e1,b2-e2,...

This function use reference to avoid copying.

The index is 1-based

=cut
sub subsequence {
    my $refseq=shift;
    my $ex=shift;
    if (!ref $ex) {
        my @arr=split /,\s*/, $ex;
        @arr = map { [split /-/] } @arr;
        $ex=\@arr;
    }
    my $subseq;
    foreach my $x (@$ex) {
        if (!defined $x) {
            die "empty exon\n";
        }
        if ($x->[0] eq '' || $x->[1] eq '') {
            die "exon range problem inside subsequence()\n";
        }
        $subseq .= substr($$refseq, $x->[0]-1, $x->[1]-$x->[0]+1);
    }
    return \$subseq;
}

=head2 splitFastaSeqfile($inputfile, $binsize, $numfilesperdir, $replace)

    parameters
      $numfilesperdir if set this number to 0, nothing will happen
      if set to >0 then it will pack the sequence files
        into directoires with $numfilesperdir files
        per directory.
      $replace [0, 1] to replace existing files or not
        default is to replace.

 @param $inputfile could be a file in the current directory
or a path to a file in a different directory.

 @param $binsize is the number of sequences per file.

 @param $numfilesperdir is an optional parameter, if given
then the files will be put into directories, with
$numfilesperdir

    if $numfilesperdir is given 0 then nothing 
    will happen.

Should produce results in the current directory
even if the input file has a path info
return the number of files

 @return a reference to an array containing all the files 
or directory names produced depending whether the 
$numfilesperdir is given or not.

The file names follow the following format depends on the 
input file format:
    1. Input file has no extension, the the output file
       will be infile#
    2. Input file has extension infile.fas or infile.fasta
       output file will be named as infile#.fas

=cut
sub splitFastaSeqfile {
    my $inputfile=shift;
    my $binsize=shift;
    my $dirsize=shift;
    my $rep=shift;
    if (! defined $rep) { $rep=1; } # default is to replace.

    if (!$binsize) { $binsize=300; } 
                  # default binsize is 300 seq per file
    open IN, "<$inputfile" or die $!, " cannot open fasta input file $inputfile\n";
    my $count=0;
    my $fext;
    my ($bin, $oufile);
    my $iname = fileOfPath($inputfile);
    my @filesProduced = ();
    if ($iname =~ /(.+)\.fas(?:ta)?$/) {
        $iname = $1;
        $fext="fas";
    }
    use integer;

    $_ = <IN>;
    while ($_ && /^>/) {
        # 0 is not allowed in SGE array job ids
        $bin = $count/$binsize + 1;
        $oufile = $iname . $bin;
        if ($fext) { $oufile .= ".$fext"; }
        push @filesProduced, $oufile;
        my $dowrite = 0;
        if (! -e $oufile || -z $oufile || $rep) {
            open OU, ">$oufile";
            $dowrite=1;
        }
        do {
            s/\t/ /g;
            if ($dowrite) {
                print OU;  # print >seqname line
            }
            $_ = <IN>;
            while ($_ && !/^>/) {  # print the sequence itself
                if ($dowrite) {
                    print OU;
                }
                $_ = <IN>;
            }
            #if (!$_) { last; }     # eof reached; stop
            $count++;
        } while ($_ && $count % $binsize != 0);
        close OU;
    }
    close IN;
    close OU;
    print STDERR "$count sequences in $inputfile split into ", $bin, " files in the current directory\n";
    #return $bin;  # the last file extension number
    my $infilestem;
    if ($inputfile =~ /(.+)\.(?:.+)/) {
        $infilestem = $1;
    }
    else {
        $infilestem = $inputfile;
    }
    my @dirsProduced;
    if ($dirsize) {
        use integer;
        my $i=0;
        while ($i < @filesProduced) {
            my $d = $i/$dirsize + 1;
            my $dirname = $infilestem . "DIR" . $d;
            mkdir $dirname;
            push @dirsProduced, $dirname;
            while ($i < @filesProduced &&
                ($i/$dirsize + 1) == $d)
            {
                system("mv", $filesProduced[$i], $dirname);
                ++$i;
            }
        }
        print STDERR scalar(@dirsProduced), " directories created\n";
        return \@dirsProduced;
    }
    else {
        return \@filesProduced;
    }
}

=head2 countFasta
 count the total number of sequences in a fasta file
=cut
sub countFasta {
   my $file=shift;
   open IN, "<$file" or die $file;
   my $cnt=0;
   while (<IN>) {
      if (/^>/) { ++$cnt; }
   }
   print STDERR "$cnt sequences in $file\n";
   close IN;
   return $cnt;
}

=head2 breakfasta
   @param($infile, $numfiles)

 break fasta file $infile into $numfiles smaller files

=cut
sub breakfasta {
   my $file=shift; # inputfile
   my $nof=shift; # number of files
   print STDERR "breaking $file into $nof pieces\n";
   #my $rep=shift;
   #my $fext=shift;
   my $numseq=countFasta($file);
   my $binsize=int($numseq/$nof) + 1;
   open IN, "<$file" or die $!;
   $_ = <IN>;
   my ($bin, $outfile, $count, $fext);
   my @filesProduced=();
   $fext='fas';
   $count=0;
   my $fstem=fileStem($file);
   while ($_ && /^>/) {
      # 0 is not allowed in SGE array job ids, 
      # needs to start from 1
      if ($count % $binsize == 0) {
         $bin = $count/$binsize + 1;
         $oufile = $fstem . $bin;
         if ($fext) { 
               $oufile .= ".$fext"; 
         }
         push @filesProduced, $oufile;
         open OU, ">$oufile" or die "could not open $oufile inside Bioseq::breakfasta!\n";
         if ($bin % 50 == 0) {
            print STDERR "working on $bin ...\n";
         }
      }
      s/\t/ /g; # replace tab with space
      print OU;  # print >seqname line
      ++$count;
# read and write sequence
      $_ = <IN>;
      while ($_ && !/^>/) {  # print the sequence itself
         print OU;
         $_ = <IN>;
      }
   }
   return \@filesProduced;
}

# fileOfPath extract the file name from 
# the full path to the file
sub fileOfPath {
    my $path = shift;
    $path =~ s/\/$//;
    my @arr = split(/\//, $path);
    if (@arr > 1) { return $arr[$#arr]; }
    return $path;
}

sub getFasPrimaryId {
# take a nr type compound id and extract the first ID
    my $longid = shift;
    chomp $longid;
    if ($longid =~ /gi\|(.+?)\|/) {
        return $1;
    }
    elsif ($longid !~ /\|/ && $longid !~ /\s/) {
        return $longid;
    }
    else {
        die "$longid contains not primary id, modify bioseq.pm!\n";
    }
}

# given a header that may contain 
# assume only one sequence per header
# the caller is responsible to process
# multiple headers
sub getFasAllId {
    my $head = shift;
    $head =~ s/^>//;
    if ($head =~ /(.+?) /) { $head = $1; }

    $head =~ s/\|$//;
    $head =~ s/\|+/|/;
    my @arr = split /\|/, $head;
    my @ids = ();
    foreach my $e (@arr) {
        if (length($e) > 3) {
            push @ids, $e;
        }
    }
    return @ids;
}

# also extract the source database
# fasta nr does not have any format
# any effort trying to get this information is wasted
sub getFasAllId_dbsrc {
    my $head = shift;
    #print STDERR "processing $head inside getFasAllId_dbsrc ...\n";
    $head =~ s/^>//;
    if ($head =~ /(.+?) /) { $head = $1; }

    $head =~ s/\|$//;
    $head =~ s/\|\|/|/;
    #print STDERR "clean head: $head\n";

    my @arr = split /\|/, $head;
    my @ids = ();
    my $i=0;
    my $db='';
    while ($i<@arr) {
        if ($arr[$i] eq 'gi' || $arr[$i] eq 'gb'
                || $arr[$i] =~ /^[a-z]+$/ ) 
        { # db name
            $db=$arr[$i++];
        }
        if (length $arr[$i] > 1) {
            push @ids, [$arr[$i], $db];
        }
        $i++;
    }

    return @ids;
}

=head2 fileStem($filename)

    return the first string before the period "." of the file name
    if the file name does not have period then return the file

For examle fileStem("protein.fas") is protein
           fileStem("somefile") is somefile

=cut

sub fileStem {
    my $file=shift;
    if ($file =~ /(.+?)\./) {
        return $1;
    }
    else {
        return $file;
    }
}


1;

__END__

