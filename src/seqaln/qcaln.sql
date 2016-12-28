/** how to filter out bad alignments
This section use the proteinSW table.
The output format of this set of program should be used.
*/

select * from proteinSWX limit 5;
-- n is a parameter
-- alnlen includes gaps, we don't want gaps
select identical/alnlen as identity, 
	alnlen, E, score,
	(1 + 480*pow(
	            (alnlen-hitins-modelins), 
				-0.32*(1 + exp( (hitins+modelins-alnlen)/1000 ) )
				))/100 as cutoff
from modelaln
where quality>=0 and identical/alnlen < 0.01*(1 + 480*pow((alnlen-hitins-modelins), -0.32*(1 + exp( (hitins+modelins-alnlen)/1000 ) ) ));

-- below the cutoff is marked as -1

update modelaln set quality=-1 
where quality>=0 and 
	identical/alnlen < 
		4.8*pow((alnlen-hitins-modelins), -0.32*(1 + exp((hitins+modelins-alnlen)/1000)));

-- quality=0 when at border line
update modelaln set quality=0 
where quality>0 and 
	identical/alnlen < 
		0.02+4.8*pow((alnlen-hitins-modelins), -0.32*(1 + exp((hitins+modelins-alnlen)/1000)));

-- compute entropy of the alignment
-- use the alignment package
