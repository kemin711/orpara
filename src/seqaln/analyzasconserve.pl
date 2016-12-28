### perl ###
##!/usr/bin/perl -w

use JGIDB;
use strict;

my ($host, $database, $user, $password);
$host='frylock-db';
$database='kzcombest';
my $crosstab='ascross';

my $db=JGIDB->new;
$db->connect(host =>  $host, database => $database, user => $user, password => $password);

createCrossTable($crosstab);
# tablename, db1, db1, db1alias, db2alias
my @input=(['Aspac1xSpoth2alnbest', 'Aspac1_1', 'Spoth2', 'Aspac1', 'Spoth2'],
   ['Aspca3xAspac1alnbest', 'Aspca3', 'Aspac1_1', 'Aspca3', 'Aspac1'],
   ['Lacbi2xAgabi2alnbest', 'Lacbi2', 'Agabi_varbisH97_2', 'Lacbi2', 'Agabi2'],
   ['Lacbi2xAspac1alnbest', 'Lacbi2', 'Aspac1_1', 'Lacbi2', 'Aspac1']);
foreach my $row (@input) {
   my $gp=bestaln2GenePair($row->[0], $row->[1], $row->[2]);
   addtoCross($gp, $row->[3], $row->[4]);
}

my $crossexpect=createASFractionTable('asfraction');
my $sqlstr=<<ENDQ;
update $crosstab c, $crossexpect e
set c.fraction_expect=e.fraction
where c.dbname1=e.dbname1 and c.dbname2=e.dbname2
    and c.hasAS1=e.hasAS1 and c.hasAS2=e.hasAS2
ENDQ
$db->doordie($sqlstr);
dumpResult();

##################################################

sub dumpResult {
   my $sql=<<ENDQ;
select dbname1, dbname2, hasAS1, hasAS2,
   observe, fraction_observe as frob,
   fraction_expect as frex,
   fraction_observe/fraction_expect as oe
from $crosstab
order by dbname1, dbname2, hasAS1, hasAS2
ENDQ

   open OU, ">ASConserve.tab" or die $!;
   $db->showQuery($sql, *OU);
   close OU;
   print "ASConserve.tab contains the alternative splicing conservation info\n";
}

sub bestaln2GenePair {
	my $bestaln=shift; # bestaln table, such as Lacbi2xAgabi2alnbest
	my $dbname1=shift;
	my $dbname2=shift;
    print "working on $bestaln ...\n";

   # 1. clean up, second round had no effect
   my $sqlstr=<<ENDQ;
delete from $bestaln
where (seq1end-seq1begin+1)/seq1len < 0.5
    and (seq2end-seq2begin+1)/seq2len < 0.5
ENDQ
   $db->doordie($sqlstr);
   $sqlstr=<<ENDQ;
delete from $bestaln
where ((seq1end-seq1begin+1)/seq1len < 0.15
      or (seq2end-seq2begin+1)/seq2len < 0.15)
   and alnlen < 70
ENDQ
   $db->doordie($sqlstr);

   # 2. add index
   $sqlstr=<<ENDQ;
alter table $bestaln
   add primary key(seq1name,seq2name),
   add key seq1(seq1name),
   add key seq2(seq2name)
ENDQ
   $db->do($sqlstr);

   # 3. build genepair table
   my $genepair=genepairName($bestaln);
   $db->doordie("drop table if exists $genepair");
   $sqlstr=<<ENDQ;
create table $genepair as
select distinct g1.geneid geneid1, g2.geneid as geneid2
from $dbname1.combest_good g1 join $bestaln a
    on g1.modelid=a.seq1name
   join $dbname2.combest_good g2
      on g2.modelid=a.seq2name
ENDQ
   $db->doordie($sqlstr);
   $db->doordie("alter table $genepair add column nmod1 integer, add column nmod2 integer");
   $sqlstr=<<ENDQ;
update $genepair p, $dbname1.combest_good_gene g
set p.nmod1=g.nummodel
where p.geneid1=g.geneid
ENDQ
   $db->doordie($sqlstr);
   $sqlstr=<<ENDQ;
update $genepair p, $dbname2.combest_good_gene g
set p.nmod2=g.nummodel
where p.geneid2=g.geneid
ENDQ
   $db->doordie($sqlstr);
   print "$bestaln converted to $genepair\n";
   return $genepair;
}

sub addtoCross {
   my $gptab=shift;
   my $alias1=shift;
   my $alias2=shift;

   my $sqlstr=<<ENDQ;
insert into $crosstab
(dbname1, dbname2, hasAS1, hasAS2, observe)
select '$alias1' as dbname1, '$alias2' as dbname2,
   nmod1>1 as hasAS1, nmod2>1 as hasAS2, count(*) as observe
from $gptab
group by hasAS1, hasAS2
ENDQ
   $db->doordie($sqlstr);
   my $total=$db->firstRow("select sum(observe) from ascross where dbname1='$alias1' and dbname2='$alias2'");
   $total=$total->[0];
   print "$total in $gptab\n";
   $sqlstr=<<ENDQ;
update $crosstab set fraction_observe=observe/$total
where dbname1='$alias1' and dbname2='$alias2'
ENDQ
   $db->doordie($sqlstr);
}

sub createCrossTable {
   my $tab=shift;
   $db->doordie("drop table if exists $tab");
# the higher the precision the higher the error!
   # don't use double!
   my $sqlstr=<<ENDQ;
create table $tab (
   dbname1 varchar(36),
   dbname2 varchar(36),
   hasAS1 boolean,
   hasAS2 boolean,
   observe integer,
   fraction_observe float,
   fraction_expect float
)
ENDQ
   $db->doordie($sqlstr);
}

sub createASFractionTable {
   my $tab=shift;

   $db->doordie("drop table if exists $tab");
   my $sqlstr=<<ENDQ;
create table $tab (
   dbname varchar(36),
   hasAS boolean,
   count integer,
   fraction float
)
ENDQ
   $db->doordie($sqlstr);
   my @input = ( 
      ['Aspac1', 'Aspac1_1'], 
      ['Agabi2', 'Agabi_varbisH97_2'],
      ['Lacbi2', 'Lacbi2'], 
      ['Aspca3', 'Aspca3'], 
      ['Spoth2', 'Spoth2']);

   foreach my $row (@input) {
      my $alias=$row->[0];
      my $dbname=$row->[1];

      $sqlstr=<<ENDQ;
insert into $tab
select '$alias' as dbname, a.hasAS, a.count, a.count/b.total as fraction
from (select nummodel > 1 as hasAS, count(*) as count
      from $dbname.combest_good_gene
      group by hasAS) a join (
         select count(*) as total
         from $dbname.combest_good_gene) b
ENDQ
      $db->doordie($sqlstr);
   }
   # create expected cross

   my $expected_cross=$tab . "_expect";
   $db->do("drop table $expected_cross");
   $db->do("drop view $expected_cross");

   $sqlstr=<<ENDQ;
create table $expected_cross (
   dbname1 varchar(36),
   dbname2 varchar(36),
   hasAs1 boolean,
   hasAS2 boolean,
   fraction float
)
ENDQ
   $db->doordie($sqlstr);
   $sqlstr=<<ENDQ;
insert into $expected_cross
select f1.dbname as dbname1, f2.dbname as dbname2,
    f1.hasAS as hasAS1, f2.hasAS as hasAS2,
        f1.fraction*f2.fraction as fraction
from $tab f1 join $tab f2
where f1.dbname != f2.dbname
ENDQ
   $db->doordie($sqlstr);
   print "$expected_cross created\n";
   return $expected_cross;
}

sub genepairName {
   my $bestaln=shift;
   if ($bestaln =~ /(.+)alnbest$/) {
      return $1 . 'genepair';
   }
   else {
      print $bestaln, " not the right format\n";
      die;
   }
}
