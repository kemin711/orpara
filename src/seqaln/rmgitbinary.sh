#/bin/sh


binary="alnexon alnglobal dbaln_fungalbt alnlocalmany analyzasconserve chimeradetect_mysql dbaln dbaln_cds dbaln_chimera dbaln_chimerarmsimp dbaln_chimerathread dbaln_fungal alnlocal dbaln_fungalrmsimp dbaln_nr exsamseq gwseedaln loadaln2mysql pickbestaln testmysql testscoremethod testaln testbraboualn"

for b in $binary; do
	git rm -f $b
done
