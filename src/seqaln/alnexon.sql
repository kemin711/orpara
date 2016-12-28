drop table alnexonAspca3Spoth1;
create table alnexonAspca3Spoth1 (
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
load data local infile 'exonalnres.tab' into table 
alnexonAspca3Spoth1 ignore 1 lines;

drop view genepair_exaln_quality;
create view genepair_exaln_quality as
select Aspca3gene, Spoth1gene, e.modid1, e.modid2, 
	a.numod1, a.numod2,
	min(e.alnquality) as quality
from bestalnAspca3Spoth1 a join alnexonAspca3Spoth1 e
	on a.seq1name=e.modid1 and a.seq2name=e.modid2
group by Aspca3gene, Spoth1gene
;

select quality, count(*)
from genepair_exaln_quality
group by quality;

+---------+----------+
| quality | count(*) |
+---------+----------+
|       0 |      335 | 
|       1 |      723 | 
|       2 |     1255 | 
|       3 |     1909 | 
|       4 |     2226 | 
+---------+----------+

-- AS set
select quality, count(*)
from genepair_exaln_quality
where numod1>1 and numod2>1
group by quality;

+---------+----------+
| quality | count(*) |
+---------+----------+
|       0 |      179 | 
|       1 |      211 | 
|       2 |      518 | 
|       3 |      930 | 
|       4 |      458 | 
+---------+----------+

-- non AS set
select numod1>1 as AS1, numod2>1 as AS2, quality, count(*) as count
from genepair_exaln_quality
group by AS1, AS2, quality;

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




