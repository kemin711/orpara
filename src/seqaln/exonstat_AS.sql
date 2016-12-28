/** exon statistics with reagards to AS 
	AS genes have shorter average exons compared to non-AS genes.
*/
select avg(exlen)
from (
	select m.exonLength/numexon as exlen
	from combest_good_gene g join combest_good m on g.repmod=m.modelid
	where g.nummodel > 1
	) foo
;
-- 656

select avg(exlen)
from (
	select m.exonLength/numexon as exlen
	from combest_good_gene g join combest_good m on g.repmod=m.modelid
	where g.nummodel = 1
	) foo
;
-- 869
-- this can be derived from correlation for AS x num exons
