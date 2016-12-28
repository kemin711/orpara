## will directory will only build binaries. 
## The uppler directory Makefile will pick the header and lib files
## and added into the library

bin_PROGRAMS = alnlocal alnglobal alnlocalmany alnglobalparam bbaln braboualnparam checkpalindrome

## test_BINARY=testkmer testChimera testbraboualn testscoremethod testmysql testaln 

LDADD = ../libkemin.la
AM_CXXFLAGS = -I..
AM_LDFLAGS = -L..

alnlocal_SOURCES = alnlocal.cpp scorematrix.h scorematrix.cpp dynalnt.h

alnlocalmany_SOURCES = alnlocalmany.cpp scorematrix.h scorematrix.cpp dynalnt.h

alnglobal_SOURCES = alnglobal.cpp scorematrix.h scorematrix.cpp dynalnt.h

alnglobalparam_SOURCES = alnglobalparam.cpp scorematrix.h scorematrix.cpp dynalnt.h

bbaln_SOURCES = bbaln.cpp braboualn.h braboualn.cpp
braboualnparam_SOURCES = braboualnparam.cpp braboualn.cpp braboualn.h

checkpalindrome_SOURCES = checkpalindrome.cpp kmerhelper.cpp kmert.h kmer.h kmer.cpp
