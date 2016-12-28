drop table alnexonAgabi2Aspca3;
create table alnexonAgabi2Aspca3 (
	modid1 integer unsigned,
	modid2 integer unsigned, 
	exlen1 text, exlen2 text, 
	exaln text, 
	numatch integer, numismatch integer, numindel integer, 
	alnbothends boolean, 
	alnquality integer,
	primary key(modid1, modid2),
	key m1id(modid1), key m2id(modid2)
);
load data local infile 'Agabi2Aspca3aln.xln.tab' into table 
alnexonAgabi2Aspca3 ignore 1 lines;

drop view if exists genepair_exaln_quality_Agabi2Aspca3;
create view genepair_exaln_quality_Agabi2Aspca3 as
select Agabi2gene, Aspca3gene, e.modid1, e.modid2, 
	a.numod1, a.numod2,
	min(e.alnquality) as quality
from Agabi2Aspca3alnbest a join alnexonAgabi2Aspca3 e
	on a.seq1name=e.modid1 and a.seq2name=e.modid2
group by Agabi2gene, Aspca3gene
;

<<<<<<< alnexonAgabiAspca3.sql
=======
select quality, count(*)
from genepair_exaln_quality_Agabi2Aspca3
group by quality;

Aspca3 x Spoth1
+---------+----------+
| quality | count(*) |
+---------+----------+
|       0 |      335 | 
|       1 |      723 | 
|       2 |     1255 | 
|       3 |     1909 | 
|       4 |     2226 | 
+---------+----------+
Agabi2 x Aspca3
+---------+----------+
| quality | count(*) |
+---------+----------+
|       0 |        8 | 
|       1 |       39 | 
|       2 |      269 | 
|       3 |     1935 | 
|       4 |     2330 | 
+---------+----------+
>>>>>>> 1.3


-- non AS set
select numod1>1 as AS1, numod2>1 as AS2, quality, count(*) as count
from genepair_exaln_quality_Agabi2Aspca3
group by AS1, AS2, quality;

Aspca3 x Spoth1
+------+------+---------+-------+
| AS1  | AS2  | quality | count |
+------+------+---------+-------+
|    0 |    0 |       0 |    42 | 
|    0 |    0 |       1 |   246 | 
|    0 |    0 |       2 |   226 | 
|    0 |    0 |       3 |   222 | 
|    0 |    0 |       4 |   736 | 
|    0 |    1 |       0 |    68 | 
|    0 |    1 |       1 |   122 | 
|    0 |    1 |       2 |   247 | 
|    0 |    1 |       3 |   419 | 
|    0 |    1 |       4 |   595 | 
|    1 |    0 |       0 |    46 | 
|    1 |    0 |       1 |   144 | 
|    1 |    0 |       2 |   264 | 
|    1 |    0 |       3 |   338 | 
|    1 |    0 |       4 |   437 | 
|    1 |    1 |       0 |   179 | 
|    1 |    1 |       1 |   211 | 
|    1 |    1 |       2 |   518 | 
|    1 |    1 |       3 |   930 | 
|    1 |    1 |       4 |   458 | 
+------+------+---------+-------+

Agabi2 X Aspca3
+------+------+---------+-------+
| AS1  | AS2  | quality | count |
+------+------+---------+-------+
|    0 |    0 |       0 |     2 | 
|    0 |    0 |       1 |    14 | 
|    0 |    0 |       2 |    66 | 
|    0 |    0 |       3 |   298 | 
|    0 |    0 |       4 |   709 | 
|    0 |    1 |       0 |     5 | 
|    0 |    1 |       1 |     9 | 
|    0 |    1 |       2 |    64 | 
|    0 |    1 |       3 |   667 | 
|    0 |    1 |       4 |   699 | 
|    1 |    0 |       1 |     4 | 
|    1 |    0 |       2 |    37 | 
|    1 |    0 |       3 |   260 | 
|    1 |    0 |       4 |   464 | 
|    1 |    1 |       0 |     1 | 
|    1 |    1 |       1 |    12 | 
|    1 |    1 |       2 |   102 | 
|    1 |    1 |       3 |   710 | 
|    1 |    1 |       4 |   458 | 
+------+------+---------+-------+

create table exaln_Agabi2Aspca3_summary as
select numod1>1 as AS1, numod2>1 as AS2, quality, count(*) as count
from genepair_exaln_quality_Agabi2Aspca3
group by AS1, AS2, quality;

-- to make a pivot table:
create view exaln_Agabi2Aspca3_summarypivot as
select distinct AS1, AS2, 
	(select count from exaln_Agabi2Aspca3_summary s 
		where s.AS1=q.AS1 and s.AS2=q.AS2 and s.quality=0) as perfect,
	(select count from exaln_Agabi2Aspca3_summary s2
		where s2.AS1=q.AS1 and s2.AS2=q.AS2 and s2.quality=1) as nearperfect,
	(select count from exaln_Agabi2Aspca3_summary s3
		where s3.AS1=q.AS1 and s3.AS2=q.AS2 and s3.quality=2) as indel,
	(select count from exaln_Agabi2Aspca3_summary s4
		where s4.AS1=q.AS1 and s4.AS2=q.AS2 and s4.quality=3) as partial,
	(select count from exaln_Agabi2Aspca3_summary s5
		where s5.AS1=q.AS1 and s5.AS2=q.AS2 and s5.quality=4) as none
from exaln_Agabi2Aspca3_summary q
;

mysql -h genome-db -e "select * from exaln_Agabi2Aspca3_summarypivot" kzcombest > exaln_AS_Agabi2Aspca3.tab

