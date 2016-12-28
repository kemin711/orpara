### perl ###
##!/usr/bin/perl -w

use DBI;
use JGIDB;

# load aln results to mysql database, simple script
my ($user,$password,$host,$database);
$host='genome-db';
$database='kzcombest';
#my $alntab="Aspca3Spoth1bestaln";
#my $alntabfile="aln.tab";
my ($alntabfile, $alntab);
my $i=0;
while ($ARGV[$i]) {
    if ($ARGV[$i] eq "-h") { $host=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq '-d') { $database=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq '-f') { $alntabfile=$ARGV[++$i]; }
    else {
        $alntab=$ARGV[$i];
    }
    ++$i;
}
if (!$alntabfile || !$alntab) {
    usage();
}

my $db=JGIDB->new;
$db->connect(user=>$user, password=>$password, host=>$host, database=>$database);
if (!$user) { $user=$db->getUser(); }
if (!$password) { $password=$db->getPassword(); }

load_table($alntab, $alntabfile);
print "$alntab loaded from $alntabfile\n";

##################################################333
sub usage {
    print STDERR "Usage: loadaln2mysql -h <host> -d <database> -f <alntabfile> <tablename>\n";
    exit(1);
}

sub load_table {
    my $tab=shift;
    my $file=shift;
    create_alntab($tab);
    my $cmd="load data local infile '$file' into table $tab ignore 1 lines";
    $db->doordie($cmd);
}

sub create_alntab {
    my $tab=shift;
    $db->doordie("drop table if exists $tab");
    my $sqstr=<<ENDQ;
create table $tab (
    seq1name char(36),
    seq1len  integer,
    seq2name char(36),
    seq2len integer,
    score integer, 
    identical integer, similar integer, alnlen integer,
    seq1numgap integer, seq2numgap integer, 
    seq1gaplen integer, seq2gaplen integer,
    seq1begin integer, seq1end integer,
    seq2begin integer, seq2end integer,
    seq1entropy float, seq2entropy float,
    primary key(seq1name, seq2name),
    key seq1id(seq1name),
    key seq2id(seq2name)
)
ENDQ
    $db->doordie($sqstr);
}
