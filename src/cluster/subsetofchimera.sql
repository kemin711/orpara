/** use nrep chimera to obtain 
chimera in particular tracks, such as the GeneCatalog track
*/

select *
from chimera_nrepmodels;

select sfStarts, sfEnds,
	trim(trailing ',' from sfStarts) as trimmed
from GeneCatalog;

create table GeneCatalog_Jan16 as select * from GeneCatalog;

update GeneCatalog 
set sfStarts=trim(trailing ',' from sfStarts),
	sfEnds=trim(trailing ',' from sfEnds)
;

select count(*)

drop table if exists chimeraGeneCatalog;
drop table if exists chimera_GeneCatalog;
create table chimera_GeneCatalog as
select c.*
from chimera_nrepmodels c join GeneCatalog g
	on c.chrom=g.chrom and c.strand=g.strand
	and c.sfstarts=g.sfStarts and
	c.sfends=g.sfEnds
;

select count(*) from chimera_GeneCatalog c join chimerasnrep_web w
	on c.name=w.name
;

select distinct colorR, colorG, colorG
from chimerasnrep_web;

update chimerasnrep_web w, chimera_GeneCatalog c
set w.colorR=254, w.colorG=1, w.colorB=1
where w.name=c.name;



drop table if exists chimera_2_GeneCatalog_mapping;
create table chimera_2_GeneCatalog_mapping as
select c.name as chimera_name,
	g.name as GC_name,
	c.proteinid as chimera_pid,
	g.proteinid as GC_pid
from chimera_nrepmodels c join GeneCatalog g
	on c.chrom=g.chrom and c.strand=g.strand
	and c.sfstarts=g.sfStarts and
	c.sfends=g.sfEnds
;

select m.*, p1.seq as chimera_pep, p2.seq as GC_pep
from chimera_2_GeneCatalog_mapping m join protein p1
	on m.chimera_pid=p1.proteinId
	join protein p2 on m.GC_pid=p2.proteinId
where m.chimera_pid != m.GC_pid
	and p1.seq != p2.seq
;


