/** input table modelaln_nrep */
-- output table modelaln_forchimera

create table modelaln_nrepmaxq as 
select id, hittaxid, max(score) as maxscore
from modelaln_nrep
where quality >=0
group by id, hittaxid;
