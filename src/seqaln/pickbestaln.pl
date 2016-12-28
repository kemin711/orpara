### perl ###
##!/usr/bin/perl -w

use DBI;
use JGIDB;
use strict;

# pick the best in both directions

my ($user,$password,$database,$host, $alntab);
$database='kzcombest';
$host='frylock-db';
my $cutoff=1; # best 90% of the top
my $i=0;
while ($ARGV[$i]) {
	if ($ARGV[$i] eq "-h") { $host=$ARGV[++$i]; }
	elsif ($ARGV[$i] eq "-d") { $database=$ARGV[++$i]; }
	elsif ($ARGV[$i] eq "-c") { $cutoff=$ARGV[++$i]; }
	elsif ($ARGV[$i] eq "?" || $ARGV[$i] eq '--help') { usage(); }
    else {
        $alntab=$ARGV[$i];
    }
    ++$i;
}
if (!$database) { usage(); }

my $db=JGIDB->new;
$db->connect(host=>$host, database=>$database, user=>$user, password=>$password);
pickbest($alntab);

##############################################

sub usage {
    print STDERR "Usage: pickbestaln -h <host> -d <database> <alntab>\n",
        " This program will produce alntabbest\n",
        " Options -c <cutoff> default 1, you can specify a value between\n",
        " (0,1]\n";
    exit(1);
}


=head2 pickbest

This can be made into a generic function using 
two columns as input

=cut
sub pickbest {
    my $tab=shift;
    my $qmax=$tab . "_qmax";
    print STDERR "creating $qmax table ...\n";
    $db->doordie("drop table if exists $qmax");
    my $sqlstr=<<ENDQ;
create table $qmax as
select seq1name, max(score) as maxscore
from $tab group by seq1name
ENDQ
    $db->doordie($sqlstr);
    $db->doordie("alter table $qmax add primary key(seq1name)");

    my $qbest=$tab . "_qbest";
    print STDERR "creating $qbest table ...\n";
    $db->doordie("drop table if exists $qbest");
    $sqlstr=<<ENDQ;
create table $qbest as
select b.*
from $tab b join $qmax q using(seq1name)
where b.score >= q.maxscore*$cutoff
ENDQ
    $db->doordie($sqlstr);

    my $tmax=$tab . "_tmax";
    print STDERR "creating $tmax table ...\n";
    $db->doordie("drop table if exists $tmax");
    $sqlstr=<<ENDQ;
create table $tmax as
select seq2name, max(score) as maxscore
from $qbest
group by seq2name
ENDQ
    $db->doordie($sqlstr);
    $db->doordie("alter table $tmax add primary key(seq2name)");

    my $best=$tab . "best";
    print STDERR "creating best blast table $best ...\n";
    $db->doordie("drop table if exists $best");
    $sqlstr=<<ENDQ;
create table $best as
select b.*
from $qbest b join $tmax m using(seq2name)
where b.score >= m.maxscore*$cutoff
ENDQ
    $db->doordie($sqlstr);
    print STDERR $best, " created\n";
}
