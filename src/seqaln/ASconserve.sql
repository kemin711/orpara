drop table if exists Aspca3Spoth1genepair;
create table Aspca3Spoth1genepair as
select distinct g1.geneid Aspca3gene, g2.geneid as Spoth1gene
from Aspca3goodshort g1 join Aspca3Spoth1alnbest a
    on g1.modelid=a.seq1name
	join Spoth1goodshort g2
	on g2.modelid=a.seq2name
;
-- 6448
select g1.geneid Aspca3gene, g2.geneid as Spoth1gene,
	max(score) as score
from Aspca3goodshort g1 join Aspca3Spoth1alnbest a
    on g1.modelid=a.seq1name
	join Spoth1goodshort g2
	on g2.modelid=a.seq2name
group by Aspca3Gene, Spoth1gene
;

-- add extra information from gene table
drop view Aspca3Spoth1genepair_nmod;
create view Aspca3Spoth1genepair_nmod as
select p.Aspca3gene, g1.nummodel as nmod1, p.Spoth1gene, g2.nummodel as nmod2
from Aspca3Spoth1genepair p join Aspca3goodgene g1
    on p.Aspca3gene=g1.geneid
	join Spoth1goodgene g2
	on p.Spoth1gene=g2.geneid
;

select nmod1, nmod2, count(*)
from Aspca3Spoth1genepair_nmod
group by nmod1, nmod2;

mysql -h genome-db -e 'select nmod1, nmod2 from Aspca3Spoth1genepair_nmod' kzcombest > Aspca3Spoth1gpAS.tab
mysql -h genome-db -e 'select nmod1, nmod2, count(*) as count from Aspca3Spoth1genepair_nmod group by nmod1,nmod2' kzcombest > Aspca3Spoth1AScount.tab



-- add more information, score, identity, coverage

create table Aspca3Spoth1ASmatch as
select nmod1 = 1 as single1, nmod2 = 1 as single2, count(*) as count
from Aspca3Spoth1genepair_nmod
group by single1, single2;

+---------+---------+-------+
| single1 | single2 | count |
+---------+---------+-------+
|       0 |       0 |  2296 | 
|       0 |       1 |  1229 | 
|       1 |       0 |  1451 | 
|       1 |       1 |  1472 | 
+---------+---------+-------+

-- this needs to match the peptide dump table, I used
-- > 20 aa cutoff!
select nummodel = 1 as single, count(*) as count,
	count(*)/(select count(*) from Agabi2goodgene_dumpped) as frac
from Agabi2goodgene_dumpped
group by single;

+--------+-------+
| single | count |
+--------+-------+
|      0 |  2966 |  AS   0.2810
|      1 |  7589 |  NoAS 0.7190
+--------+-------+

select nummodel = 1 as single, count(*) as count,
	count(*)/(select count(*) from Aspca3goodgene) as frac
from Aspca3goodgene
group by single;

+---------+----------+
| singlex | count(*) |
+---------+----------+
|       0 |     5085 | AS    0.3344
|       1 |    10122 | NoAS  0.6656
+---------+----------+
-- 33.44% with Alternative

-- Spoth1
select nummodel = 1 as single, count(*) as count,
	count(*)/(select count(*) from Spoth1goodgene) as frac
from Spoth1goodgene
group by single;

+--------+-------+--------+
| single | count | frac   |
+--------+-------+--------+
|      0 |  4115 | 0.3286 | 
|      1 |  8408 | 0.6714 | 
+--------+-------+--------+


create table Agabi2AS (
	altsp varchar(8),
	frac float );

insert into Agabi2AS values ('AS', 0.2810),('NoAs', 0.7190);

create table Aspca3AS (
	altsp varchar(8),
	frac float );

insert into Aspca3AS values ('AS', 0.3344),('NoAs', 0.6656);

create table Spoth1AS (
	altsp varchar(8),
	frac float );

insert into Spoth1AS values ('AS', 0.3286),('NoAs', 0.6714);

create table Aspca3Spoth1AStest as
select concat(a1.altsp, ' x ', a2.altsp) as ASmatch,
	a1.frac*a2.frac as expected
from Aspca3AS a1 join Spoth1AS a2;

alter table Aspca3Spoth1AStest add column observe integer;
-- use Aspca3Spoth1Match
select * from Aspca3Spoth1ASmatch;
+---------+---------+-------+
| single1 | single2 | count |
+---------+---------+-------+
|       0 |       0 |  2296 | 
|       0 |       1 |  1229 | 
|       1 |       0 |  1451 | 
|       1 |       1 |  1472 | 
+---------+---------+-------+

update Aspca3Spoth1AStest set observe=2296
where ASmatch='AS x AS';
update Aspca3Spoth1AStest set observe=1229
where ASmatch='AS x NoAS';
update Aspca3Spoth1AStest set observe=1451
where ASmatch='NoAS x AS';
update Aspca3Spoth1AStest set observe=1472
where ASmatch='NoAS x NoAS';

+---------+---------+-------+
| single1 | single2 | count |
+---------+---------+-------+
|       0 |       0 |  1283 | 
|       0 |       1 |   765 | 
|       1 |       0 |  1444 | 
|       1 |       1 |  1089 | 
+---------+---------+-------+

select * from Aspca3Spoth1AStest;
alter table Aspca3Spoth1AStest add column expected_count integer;

select expected*(select sum(observe) from Aspca3Spoth1AStest) as expected_count
from Aspca3Spoth1AStest;

select sum(observe) from Aspca3Spoth1AStest;
-- 6448

update  Aspca3Spoth1AStest
set expected_count=expected*6448
;

select ASmatch, round(expected, 3) as expected, 
	observe as obsvCnt, expected_count as expCnt,
	observe/(select sum(observe) from Aspca3Spoth1AStest) as obsfreq
from Aspca3Spoth1AStest;
+-------------+----------+---------+--------+---------+
| ASmatch     | expected | obsvCnt | expCnt | obsfreq |
+-------------+----------+---------+--------+---------+
| AS x AS     |    0.110 |    2296 |    709 |  0.3561 | 
| NoAs x AS   |    0.219 |    1451 |   1410 |  0.2250 | 
| AS x NoAs   |    0.225 |    1229 |   1448 |  0.1906 | 
| NoAs x NoAs |    0.447 |    1472 |   2882 |  0.2283 | 
+-------------+----------+---------+--------+---------+

--- more clean measurement
alter table bestalnAspca3Spoth1 add column Aspca3gene integer unsigned, 
	add column Spoth1gene integer unsigned
;

update bestalnAspca3Spoth1 a, Aspca3goodshort g
set a.Aspca3gene=g.geneid
where a.seq1name=g.modelid;
-- 10914
update bestalnAspca3Spoth1 a, Spoth1goodshort g
set a.Spoth1gene=g.geneid
where a.seq2name=g.modelid;

optimize table bestalnAspca3Spoth1;


select Aspca3gene, Spoth1gene, seq1name, seq2name, score, identical/alnlen as identity
from bestalnAspca3Spoth1
order by Aspca3gene, Spoth1gene;

select min(identical/alnlen) from bestalnAspca3Spoth1;
-- 0.1936

alter table bestalnAspca3Spoth1 add column numod1 integer, add column numod2 integer;

alter table bestalnAspca3Spoth1 add key (Aspca3gene), add key (Spoth1gene);

update bestalnAspca3Spoth1 a, Aspca3goodgene g
set a.numod1=g.nummodel
where a.Aspca3gene=g.geneid;

update bestalnAspca3Spoth1 a, Spoth1goodgene g
set a.numod2=g.nummodel
where a.Spoth1gene=g.geneid;

select Aspca3gene, Spoth1gene, seq1name, seq2name, 
	score, identical/alnlen as identity, seq1gaplen as gap1, 
	seq1gaplen/seq1len as gap1frac, seq2gaplen as gap2, seq2gaplen/seq2len as gap2frac, 
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2
from bestalnAspca3Spoth1
where numod1>1 and numod2>1 and identical/alnlen>0.39
	and seq1gaplen < 15 and seq2gaplen < 15 
	and (seq1end-seq1begin+1)/seq1len > 0.9
	and (seq2end-seq2begin+1)/seq2len > 0.9
;


select Aspca3gene, Spoth1gene, count(*) as count
from (
select Aspca3gene, Spoth1gene, seq1name, seq2name, 
	score, identical/alnlen as identity, seq1gaplen as gap1, 
	seq1gaplen/seq1len as gap1frac, seq2gaplen as gap2, seq2gaplen/seq2len as gap2frac, 
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2
from bestalnAspca3Spoth1
where identical/alnlen>0.39
	and seq1gaplen < 15 and seq2gaplen < 20 
	and (seq1end-seq1begin+1)/seq1len > 0.9
	and (seq2end-seq2begin+1)/seq2len > 0.9
) foo 
group by Aspca3gene, Spoth1gene
;
numod1>1 and numod2>1 and 
	and numexon1=numexon2
having count > 1
--319 gene with >39% identity, > 90% coverage of both proteins
--              total gap length < 20 nt
-- identical number of exons for both matched transcript
-- without identical number of exons
-- 732
-- outof total of genes showing any match without 
-- condition of numod1>1 and numod2 > 1
-- 1376 satisfy this condition
drop view Aspca3Spoth1SameAS;
create view Aspca3Spoth1SameAS as
select Aspca3gene, Spoth1gene, seq1name, seq2name, 
	score, identical/alnlen as identity, seq1gaplen as gap1, 
	seq1gaplen/seq1len as gap1frac, seq2gaplen as gap2, seq2gaplen/seq2len as gap2frac, 
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2
from bestalnAspca3Spoth1
where identical/alnlen>0.39 and numod1>1 and numod2>1 and abs(numexon1-numexon2)<2
	and seq1gaplen < 20 and seq2gaplen < 20 
	and (seq1end-seq1begin+1)/seq1len > 0.9
	and (seq2end-seq2begin+1)/seq2len > 0.9
;

select seq1name, concat(g1.genomicid, ':', g1.begin, '-', g1.end) as Aspca3loc,
	seq2name, concat(g2.genomicid, ':', g2.begin, '-', g2.end) as Spoth1loc, identity
from Aspca3Spoth1SameAS a join Aspca3goodgene g1 on a.Aspca3gene=g1.geneid
	join Spoth1goodgene g2 on a.Spoth1gene=g2.geneid
;

/** 1200 models meet the above condition, 6/20 models have identical 
 intron structure. Estimated 360 genes with identical AS
*/

select seq1name, seq2name, concat(m1.exons, ' ', m1.gcdsb, ' ', m1.gcdse) as arg1, 
	concat(m2.exons, ' ', m2.gcdsb, ' ', m2.gcdse) as arg2
from Aspca3Spoth1SameAS a join Aspca3goodshort m1 
	on a.seq1name=m1.modelid
	join Spoth1goodshort m2 on a.seq2name=m2.modelid
;


/** more relaxed for computer program to look at alnexon program
*/
create view goodalnAspca3Spoth1 as
select a.seq1name, a.seq2name, m1.exons as exons1, m1.gcdsb as gcdsb1, 
	m1.gcdse as gcdse1, m2.exons as exons2, m2.gcdsb as gcdsb2, 
	m2.gcdse as gcdse2, a.Aspca3gene, a.Spoth1gene, 
	numod1, numod2, numexon1, numexon2
from bestalnAspca3Spoth1 a join Aspca3goodshort m1 on a.seq1name=m1.modelid 
	 join Spoth1goodshort m2 on a.seq2name=m2.modelid
where identical/alnlen>0.33 
	and (seq1end-seq1begin+1)/seq1len > 0.8
	and (seq2end-seq2begin+1)/seq2len > 0.8
;
select count(distinct Aspca3gene, Spoth1gene)
from goodalnAspca3Spoth1;
-- 4106 gene pairs, 
select count(*)
from goodalnAspca3Spoth1; -- 6877 alignments

/** wiht AS */
select count(distinct Aspca3gene, Spoth1gene), count(*)
from goodalnAspca3Spoth1 where numod1>1 and numod2>1;
-- 1714 gene pairs, 4010 aln pairs

/** without AS */
select count(distinct Aspca3gene, Spoth1gene), count(*)
from goodalnAspca3Spoth1 where numod1=1 or numod2=1;
-- 2392 gene pairs, 2867 aln pairs, 


/** add number of exons to combest_good_short to add one more parameter for
 * filtering. The query needs to be executed on individual databases. **/
drop table if exists combest_good_short;
create table combest_good_short (
	modelid integer unsigned primary key,
	geneid integer unsigned,
	numexon integer,
	mRNAlen integer,
	exons longtext,
	gcdsb integer,
	gcdse integer
);
insert into combest_good_short
select modelid, geneid, numexon, exonLength,
	exons, gcdsb, gcdse
from combest_good;

insert into combest_good_short 
select modelid, geneid, numexon, exonLength,
	exons, gcdsb, gcdse
from combest_model where modelid=30371;

alter table combest_good_short add key geneid(geneid);

cpytable -ih frylock-db -oh genome-db -id Aspca3 -od kzcombest combest_good_short

cpytable -ih frylock-db -oh genome-db -id Agabi_varbisH97_2 -od kzcombest combest_good_short

drop table Aspca3goodshort;
alter table combest_good_short
rename to Aspca3goodshort;

drop table Agabi2goodshort;
alter table combest_good_short
rename to Agabi2goodshort;

-- for some reason combest_good from Aspca3 is missing modelid 30371, protein
-- length 428. Could not figure out why, so I will just copy over this missing
-- entry

-- for Spoth1 we can copy it directly
drop table if exists Spoth1goodshort;
create table Spoth1goodshort (
	modelid integer unsigned primary key,
	geneid integer unsigned,
	numexon integer,
	mRNAlen integer,
	exons longtext, gcdsb integer, gcdse integer
);
insert into Spoth1goodshort
select modelid, geneid, numexon, exonLength, exons,
	gcdsb, gcdse
from Spoth1_Clone.combest_good;

--- add extra informatio to aln table for convinence
alter table bestalnAspca3Spoth1 add column numexon1 integer, 
	add column numexon2 integer;

update bestalnAspca3Spoth1 a, Aspca3goodshort c
set a.numexon1=c.numexon
where a.seq1name=c.modelid;

update bestalnAspca3Spoth1 a, Spoth1goodshort c
set a.numexon2=c.numexon
where a.seq2name=c.modelid;

