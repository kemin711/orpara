-- after running the program addentropy,
-- need to do a few post processing, because I have not
-- figured out how to enter NULL throught the C++ interface

update modelaln_nrep set scoreZ=NULL, idenZ=NULL 
where scoreZ=-999 and idenZ=-999;
optimize table modelaln_nrep;
select identical/alnlen as identity,
from modelaln_nrep where idenZ > 99999999 or scoreZ>999999999;


-- short simple sequences
select alnlen, identical/alnlen as identity 
from modelaln_nrep 
where (scoreZ<2 or idenZ<2) and quality>=0;

update modelaln_nrep set quality=-99 
where (scoreZ<2 or idenZ<2) and quality>=0;

select quality from modelaln_nrep 
where (scoreZ<2 or idenZ<2) and quality<0;

select * from modelaln_nrep 
where quality >= 0 and (qentropy<1.5 or tentropy < 1.5)
	and alnlen<70 and identical/alnlen<0.8;

update modelaln_nrep set quality=-95
where quality >= 0 and (qentropy<1.5 or tentropy < 1.5)
	and alnlen<70 and identical/alnlen<0.8;
