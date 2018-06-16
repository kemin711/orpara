## will directory will only build binaries. 
## The uppler directory Makefile will pick the header and lib files
## and added into the library

bin_PROGRAMS = alnlocal alnglobal alnlocalmany alnglobalparam 

# useless program: bbaln braboualnparam 

## test_BINARY=testkmer testChimera testbraboualn testscoremethod testmysql testaln 

LDADD = ../libkemin.la
AM_CXXFLAGS = -I..
AM_LDFLAGS = -L..
AM_CPPFLAGS = -D MATRIXDIR="\"$(pkgdatadir)/matrix\""

alnlocal_SOURCES = alnlocal.cpp scorematrix.h scorematrix.cpp dynalnt.h

alnlocalmany_SOURCES = alnlocalmany.cpp scorematrix.h scorematrix.cpp dynalnt.h

alnglobal_SOURCES = alnglobal.cpp scorematrix.h scorematrix.cpp dynalnt.h

alnglobalparam_SOURCES = alnglobalparam.cpp scorematrix.h scorematrix.cpp dynalnt.h

#testrepeatrun_SOURCES = testrepeatrun.cpp

#bbaln_SOURCES = bbaln.cpp braboualn.h braboualn.cpp
#braboualnparam_SOURCES = braboualnparam.cpp braboualn.cpp braboualn.h

nobase_dist_pkgdata_DATA = matrix/AAcode.txt  \
matrix/BLOSUM50 matrix/BLOSUM70.50 matrix/BLOSUM100 matrix/BLOSUM50.50 \
matrix/BLOSUM75 matrix/BLOSUM100.50 matrix/BLOSUM55 matrix/BLOSUM75.50 \
matrix/BLOSUM30 matrix/BLOSUM55.50 matrix/BLOSUM80 matrix/BLOSUM30.50 \
matrix/BLOSUM60 matrix/BLOSUM80.50 matrix/BLOSUM35 matrix/BLOSUM60.50 \
matrix/BLOSUM85 matrix/BLOSUM35.50 matrix/BLOSUM62 matrix/BLOSUM85.50 \
matrix/BLOSUM40 matrix/BLOSUM62.50 matrix/BLOSUM90 matrix/BLOSUM40.50 \
matrix/BLOSUM65 matrix/BLOSUM90.50 matrix/BLOSUM45 matrix/BLOSUM65.50 \
matrix/BLOSUMN matrix/BLOSUM45.50 matrix/BLOSUM70 matrix/BLOSUMN.50 \
matrix/DAYHOFF matrix/DNA matrix/DNAelem matrix/GONNET matrix/IDENTITY \
matrix/MATCH matrix/NUC.4.4.N matrix/NUC.4.4 matrix/NUC.4.2 \
matrix/PAM110 matrix/PAM20 matrix/PAM280 matrix/PAM380 matrix/PAM470 \
matrix/PAM120 matrix/PAM200 matrix/PAM290 matrix/PAM390 matrix/PAM480 \
matrix/PAM30  matrix/PAM40 matrix/PAM490 matrix/PAM130 matrix/PAM210 \
matrix/PAM300 matrix/PAM400 matrix/PAM50 matrix/PAM140 matrix/PAM220 \
matrix/PAM310 matrix/PAM500 matrix/PAM150 matrix/PAM230 matrix/PAM320 \
matrix/PAM410 matrix/PAM60 matrix/PAM160 matrix/PAM240 matrix/PAM330 \
matrix/PAM420 matrix/PAM70 matrix/PAM250 matrix/PAM340 matrix/PAM430 \
matrix/PAM80 matrix/PAM170 matrix/PAM350 matrix/PAM440 matrix/PAM10 \
matrix/PAM180 matrix/PAM260 matrix/PAM360 matrix/PAM450 matrix/PAM90 \
matrix/PAM100 matrix/PAM190 matrix/PAM270 matrix/PAM370 matrix/PAM460 \
matrix/PAM120.cdi matrix/PAM40.cdi matrix/PAM160.cdi matrix/PAM80.cdi \
matrix/PAM250.cdi matrix/PAM200.cdi matrix/README

SUBDIRS = . test
