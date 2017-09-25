## Process this file with automake to produce Makefile.in
## the dependency order of the directories:
## libkemin.la is the corelib: util range bioseq seqaln genericTree 
##  seqana.la and tools : seqana 
##  only application tools: seqconv
## this file is used to generate the library for distribution

bin_PROGRAMS = GrishinDistance
GrishinDistance_SOURCES=GrishinDistance.cpp evolution/evolution.h

AUTOMAKE_OPTIONS = subdir-objects
only_headers = stddev.h seqaln/dynalnt.h seqaln/kmert.h seqaln/kmerbase.h seqaln/kmerset.h gentree.h evolution/evolution.h quicksortorder.h
header_forlibs = strformat.h range.h interval.h bioseq.h codon.h Model.h FPGModel.h GModel.h fastq.h seqaln/scorematrix.h seqaln/alninfo.h seqaln/alnexaminer.h seqaln/kmer.h seqaln/kmerhelper.h ncbitaxonomy.h

cppfiles = strformat.cpp range.cpp interval.cpp codon.cpp bioseq.cpp Model.cpp FPGModel.cpp fastq.cpp seqaln/scorematrix.cpp seqaln/alninfo.cpp seqaln/alnexaminer.cpp seqaln/kmer.cpp seqaln/kmerhelper.cpp ncbitaxonomy.cpp

pkginclude_HEADERS = $(only_headers) $(header_forlibs)

lib_LTLIBRARIES = libkemin.la
libkemin_la_SOURCES = $(header_forlibs) $(cppfiles)
libkemin_la_LDFLAGS = -version-info 0:0:0

SUBDIRS = . seqaln test
