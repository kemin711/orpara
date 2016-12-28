/** use string substitutation to edit this script
*/
delete from Lacbi2xAgabi2alnbest
where (seq1end-seq1begin+1)/seq1len < 0.5
    and (seq2end-seq2begin+1)/seq2len < 0.5
;

-- short alignments, observe
select seq1name, seq2name, identical/alnlen as identity,
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2,
	alnlen
from Lacbi2xAgabi2alnbest 
where (seq1end-seq1begin+1)/seq1len< 0.1
	or (seq2end-seq2begin+1)/seq2len < 0.1
	or identical/alnlen < 0.2
	or alnlen < 50
;

delete from Lacbi2xAgabi2alnbest
where ((seq1end-seq1begin+1)/seq1len < 0.15 
	or (seq2end-seq2begin+1)/seq2len < 0.15) 
	and alnlen < 70
;

alter table Lacbi2xAgabi2alnbest
add primary key(seq1name,seq2name),
add key seq1(seq1name),
add key seq2(seq2name)
;

select avg(cov1), avg(cov2), avg(identity),
    min(cov1), min(cov2), min(identity)
from (
	select seq1name, seq2name, identical/alnlen as identity,
		(seq1end-seq1begin+1)/seq1len as cov1,
			(seq2end-seq2begin+1)/seq2len as cov2
	from Lacbi2xAgabi2alnbest
) foo;

-- method A.

-- bring in the gene table to match genes
-- only if on different server.
-- we need combest_good
-- copy over combest_good (modelid, geneid) two columns only
-- go to Agabi_varbisH97_2, make combest_good_short
-- cannot copy view!
drop view combest_good_short;
create table combest_good_short as 
select modelid, geneid from combest_good;


cpytable -ih frylock-db -oh genome-db -id Agabi_varbisH97_2 -od kzcombest combest_good_gene combest_good_short
# this program cannot copy view!
# also need to copy over the peptide table use to do the blat
# if we have done any filtering. In this case I used 20 aa cutoff for protein
# to be used in blast.


select g.*
from combest_good_gene g left outer join goodpeptide p
	on g.repmod=p.modelid
where p.modelid is null;

create table combest_good_gene_dumpped as
select g.*
from combest_good_gene g join goodpeptide p
	on g.repmod=p.modelid
;

cpytable -ih frylock-db -oh genome-db -id Agabi_varbisH97_2 -od kzcombest combest_good_gene_dumpped

alter table combest_good_gene rename to Agabi2goodgene;
alter table combest_good_short rename to Agabi2goodshort;
alter table Agabi2goodshort add primary key(modelid);
alter table combest_good_gene_dumpped rename to Agabi2goodgene_dumpped;

-- method B
-- if on the same server, then can use different databases.
drop table if exists Lacbi2xAgabi2genepair;
create table Lacbi2xAgabi2genepair as
select distinct g1.geneid geneid1, g2.geneid as geneid2
from Lacbi2.combest_good g1 join Lacbi2xAgabi2alnbest a
    on g1.modelid=a.seq1name
	join Agabi_varbisH97_2.combest_good g2
	on g2.modelid=a.seq2name
;

-- 4581 pairs
alter table Lacbi2xAgabi2genepair add column nmod1 integer,
	add column nmod2 integer;

update Lacbi2xAgabi2genepair p, Lacbi2.combest_good_gene g
set p.nmod1=g.nummodel
where p.geneid1=g.geneid;

update Lacbi2xAgabi2genepair p, Agabi_varbisH97_2.combest_good_gene g
set p.nmod2=g.nummodel
where p.geneid2=g.geneid;

-- through view is not so good

create table Lacbi2xAgabi2ASmatch as
select nmod1 > 1 as hasAS1, nmod2 > 1 as hasAS2, count(*) as count
from Lacbi2xAgabi2genepair
group by hasAS1, hasAS2;

create table ascross (
	dbname1 varchar(36),
	dbname2 varchar(36),
	hasAS1 boolean,
	hasAS2 boolean,
	observe integer,
	fraction_observe float,
	fraction_expect float
);

select nmod1>1 as hasAS1, nmod2>1 as hasAS2, count(*) as count
from Lacbi2xAgabi2genepair
group by hasAS1, hasAS2;

select dbname1, dbname2, hasAS1, hasAS2,
	observe/(select sum(observe) from ascross
	where  dbname1='Lacbi2' and dbname2='Agabi2') as fraction
from ascross
where dbname1='Lacbi2' and dbname2='Agabi2'
;


insert into ascross 
(dbname1, dbname2, hasAS1, hasAS2, fraction_observe)
select 'Lacbi2' as dbname1, 'Agabi2' as dbname2, 
	hasAS1, hasAS2, 
	count/(select sum(count) from Lacbi2xAgabi2ASmatch) as fraction_observe
from Lacbi2xAgabi2ASmatch;

select dbname1, dbname2, hasAS1, hasAS2,
	observe, fraction_observe as frob,
	fraction_expect as frex, 
	fraction_observe/fraction_expect as oe
from ascross
order by dbname1, dbname2, hasAS1, hasAS2;
+---------+---------+--------+--------+------------------+-----------------+
| dbname1 | dbname2 | hasAS1 | hasAS2 | fraction_observe | fraction_expect |
+---------+---------+--------+--------+------------------+-----------------+
| Aspca3  | Aspac1  |      0 |      0 |         0.177975 |        0.381341 | 
| Aspca3  | Aspac1  |      0 |      1 |         0.360738 |         0.30787 | 
| Aspca3  | Aspac1  |      1 |      0 |        0.0618942 |         0.17196 | 
| Aspca3  | Aspac1  |      1 |      1 |         0.399393 |        0.138829 | 
| Aspac1  | Spoth2  |      0 |      0 |         0.142301 |        0.425247 | 
| Aspac1  | Spoth2  |      0 |      1 |        0.0873582 |        0.128054 | 
| Aspac1  | Spoth2  |      1 |      0 |         0.338412 |        0.343317 | 
| Aspac1  | Spoth2  |      1 |      1 |         0.431929 |        0.103383 | 
| Lacbi2  | Aspac1  |      0 |      0 |        0.0422422 |        0.269707 | 
| Lacbi2  | Aspac1  |      0 |      1 |         0.147147 |        0.217744 | 
| Lacbi2  | Aspac1  |      1 |      0 |         0.147147 |        0.283594 | 
| Lacbi2  | Aspac1  |      1 |      1 |         0.663463 |        0.228956 | 
| Lacbi2  | Agabi2  |      0 |      0 |         0.137396 |        0.346418 | 
| Lacbi2  | Agabi2  |      0 |      1 |        0.0599575 |        0.141032 | 
| Lacbi2  | Agabi2  |      1 |      0 |         0.495017 |        0.364255 | 
| Lacbi2  | Agabi2  |      1 |      1 |         0.307629 |        0.148294 | 
+---------+---------+--------+--------+------------------+-----------------+


-- use ASfraction.sql to create asfraction table
-- that has dbname, hasAS, count, fraction
create view asfraction_expect as
select f1.dbname as dbname1, f2.dbname as dbname2,
	f1.hasAS as hasAS1, f2.hasAS as hasAS2, 
	f1.fraction*f2.fraction as fraction
from asfraction f1 join asfraction f2
;
-- where f1.dbname='Lacbi2' and f2.dbname='Agabi2'
select * from asfraction_expect
order by dbname1, dbname2, hasAS1, hasAS2;

update ascross c, asfraction_expect e
set c.fraction_expect=e.fraction
where c.dbname1=e.dbname1 and c.dbname2=e.dbname2
	and c.hasAS1=e.hasAS1 and c.hasAS2=e.hasAS2
;


select * from ascross;


+---------+---------+-------+
| single1 | single2 | count |
+---------+---------+-------+
|       0 |       0 |  1283 | 
|       0 |       1 |   765 | 
|       1 |       0 |  1444 | 
|       1 |       1 |  1089 | 
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
	count(*)/(select count(*) from Lacbi2goodgene) as frac
from Lacbi2goodgene
group by single;

+---------+----------+
| singlex | count(*) |
+---------+----------+
|       0 |     5085 | AS    0.3344
|       1 |    10122 | NoAS  0.6656
+---------+----------+
-- 33.44% with Alternative


create table Agabi2AS (
	altsp varchar(8),
	frac float );

insert into Agabi2AS values ('AS', 0.2810),('NoAs', 0.7190);

create table Lacbi2AS (
	altsp varchar(8),
	frac float );

insert into Lacbi2AS values ('AS', 0.3344),('NoAs', 0.6656);

create table Lacbi2xAgabi2AStest as
select concat(a1.altsp, ' x ', a2.altsp) as ASmatch,
	a1.frac*a2.frac as expected
from Agabi2AS a1 join Lacbi2AS a2;

alter table Lacbi2xAgabi2AStest add column observe integer;
update Lacbi2xAgabi2AStest set observe=1283
where ASmatch='AS x AS';
update Lacbi2xAgabi2AStest set observe=765
where ASmatch='AS x NoAS';
update Lacbi2xAgabi2AStest set observe=1444
where ASmatch='NoAS x AS';
update Lacbi2xAgabi2AStest set observe=1089
where ASmatch='NoAS x NoAS';

+---------+---------+-------+
| single1 | single2 | count |
+---------+---------+-------+
|       0 |       0 |  1283 | 
|       0 |       1 |   765 | 
|       1 |       0 |  1444 | 
|       1 |       1 |  1089 | 
+---------+---------+-------+

select * from Lacbi2xAgabi2AStest;
alter table Lacbi2xAgabi2AStest add column expected_count integer;

select expected*(select sum(observe) from Lacbi2xAgabi2AStest) as expected_count
from Lacbi2xAgabi2AStest;

select sum(observe) from Lacbi2xAgabi2AStest;
-- 4581

update  Lacbi2xAgabi2AStest
set expected_count=expected*4581
;

-- this will not work
update  Lacbi2xAgabi2AStest
set expected_count=expected*(select sum(observe) from Lacbi2xAgabi2AStest)
;
mysql:kzcombest:genome-db>select * from Lacbi2xAgabi2AStest;
+-------------+--------------------+---------+----------------+
| ASmatch     | expected           | observe | expected_count |
+-------------+--------------------+---------+----------------+
| AS x AS     | 0.0939663955842018 |    1283 |            430 | 
| NoAs x AS   |   0.24043359263792 |    1444 |           1101 | 
| AS x NoAs   |  0.187033592733288 |     765 |            857 | 
| NoAs x NoAs |  0.478566389242268 |    1089 |           2192 | 
+-------------+--------------------+---------+----------------+

