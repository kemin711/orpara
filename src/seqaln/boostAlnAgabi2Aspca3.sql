--- more clean measurement
alter table Agabi2Aspca3alnbest add column Agabi2gene integer unsigned, 
	add column Aspca3gene integer unsigned
;

update Agabi2Aspca3alnbest a, Agabi2goodshort g
set a.Agabi2gene=g.geneid
where a.seq1name=g.modelid;
-- 10914
update Agabi2Aspca3alnbest a, Aspca3goodshort g
set a.Aspca3gene=g.geneid
where a.seq2name=g.modelid;

optimize table Agabi2Aspca3alnbest;


select Agabi2gene, Aspca3gene, seq1name, seq2name, score, identical/alnlen as identity
from Agabi2Aspca3alnbest
order by Agabi2gene, Aspca3gene;

select min(identical/alnlen) from Agabi2Aspca3alnbest;
-- 0.1936

alter table Agabi2Aspca3alnbest add key (Agabi2gene), add key (Aspca3gene);
alter table Agabi2Aspca3alnbest add column numod1 integer, add column numod2 integer;

update Agabi2Aspca3alnbest a, Agabi2goodgene g
set a.numod1=g.nummodel
where a.Agabi2gene=g.geneid;

update Agabi2Aspca3alnbest a, Aspca3goodgene g
set a.numod2=g.nummodel
where a.Aspca3gene=g.geneid;
-- take a look at the result.
select Agabi2gene, Aspca3gene, seq1name, seq2name, numod1, numod2
from Agabi2Aspca3alnbest;


select Agabi2gene, Aspca3gene, seq1name, seq2name, 
	score, identical/alnlen as identity, seq1gaplen as gap1, 
	seq1gaplen/seq1len as gap1frac, seq2gaplen as gap2, seq2gaplen/seq2len as gap2frac, 
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2
from Agabi2Aspca3alnbest
where numod1>1 and numod2>1 and identical/alnlen>0.39
	and seq1gaplen < 15 and seq2gaplen < 15 
	and (seq1end-seq1begin+1)/seq1len > 0.9
	and (seq2end-seq2begin+1)/seq2len > 0.9
;

select Agabi2gene, Aspca3gene, count(*) as count
from (
select Agabi2gene, Aspca3gene, seq1name, seq2name, 
	score, identical/alnlen as identity, seq1gaplen as gap1, 
	seq1gaplen/seq1len as gap1frac, seq2gaplen as gap2, seq2gaplen/seq2len as gap2frac, 
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2
from Agabi2Aspca3alnbest
where identical/alnlen>0.39
	and seq1gaplen < 15 and seq2gaplen < 20 
	and (seq1end-seq1begin+1)/seq1len > 0.9
	and (seq2end-seq2begin+1)/seq2len > 0.9
) foo 
group by Agabi2gene, Aspca3gene
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
select Agabi2gene, Aspca3gene, seq1name, seq2name, 
	score, identical/alnlen as identity, seq1gaplen as gap1, 
	seq1gaplen/seq1len as gap1frac, seq2gaplen as gap2, seq2gaplen/seq2len as gap2frac, 
	(seq1end-seq1begin+1)/seq1len as cov1,
	(seq2end-seq2begin+1)/seq2len as cov2
from Agabi2Aspca3alnbest
where identical/alnlen>0.39 and numod1>1 and numod2>1 and abs(numexon1-numexon2)<2
	and seq1gaplen < 20 and seq2gaplen < 20 
	and (seq1end-seq1begin+1)/seq1len > 0.9
	and (seq2end-seq2begin+1)/seq2len > 0.9
;

select seq1name, concat(g1.genomicid, ':', g1.begin, '-', g1.end) as Aspca3loc,
	seq2name, concat(g2.genomicid, ':', g2.begin, '-', g2.end) as Spoth1loc, identity
from Aspca3Spoth1SameAS a join Aspca3goodgene g1 on a.Agabi2gene=g1.geneid
	join Spoth1goodgene g2 on a.Aspca3gene=g2.geneid
;

/** 1200 models meet the above condition, 6/20 models have identical 
 intron structure. Estimated 360 genes with identical AS
*/

select seq1name, seq2name, concat(m1.exons, ' ', m1.gcdsb, ' ', m1.gcdse) as arg1, 
	concat(m2.exons, ' ', m2.gcdsb, ' ', m2.gcdse) as arg2
from Aspca3Spoth1SameAS a join Agabi2goodshort m1 
	on a.seq1name=m1.modelid
	join Aspca3goodshort m2 on a.seq2name=m2.modelid
;

select * from Agabi2Aspca3alnbest;


