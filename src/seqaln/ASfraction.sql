drop table asfraction;

create table asfraction (
	dbname varchar(36),
	hasAS boolean,
	count integer,
	fraction double
);

insert into asfraction
select 'Aspca3' as dbname, a.hasAS, a.count, a.count/b.total as fraction
from (select nummodel > 1 as hasAS, count(*) as count
from Aspca3.combest_good_gene
group by hasAS) a join (select count(*) as total from Aspca3.combest_good_gene)
 b
;

insert into asfraction
select 'Aspac1' as dbname, a.hasAS, a.count, a.count/b.total as fraction
from (select nummodel > 1 as hasAS, count(*) as count
      from Aspac1_1.combest_good_gene
group by hasAS) a join (
	select count(*) as total 
	from Aspac1_1.combest_good_gene) b
;

insert into asfraction
select 'Agabi2' as dbname, a.hasAS, a.count, a.count/b.total as fraction
from (select nummodel > 1 as hasAS, count(*) as count
      from Agabi_varbisH97_2.combest_good_gene
group by hasAS) a join (
	select count(*) as total 
	from Agabi_varbisH97_2.combest_good_gene) b
;

insert into asfraction
select 'Lacbi2' as dbname, a.hasAS, a.count, a.count/b.total as fraction
from (select nummodel > 1 as hasAS, count(*) as count
      from Lacbi2.combest_good_gene
group by hasAS) a join (
	select count(*) as total 
	from Lacbi2.combest_good_gene) b
;

insert into asfraction
select 'Spoth2' as dbname, a.hasAS, a.count, a.count/b.total as fraction
from (select nummodel > 1 as hasAS, count(*) as count
      from Spoth2.combest_good_gene
group by hasAS) a join (
	select count(*) as total 
	from Spoth2.combest_good_gene) b
;

