TESTS = testAll 
check_PROGRAMS = testAll
testAll_SOURCES = main_testall.cpp bioseq_test.cpp interval_test.cpp quicksortorder_test.cpp \
						fastq_test.cpp kmean_test.cpp randombase_test.cpp
testAll_CXXFLAGS = -I.. -I../seqaln $(gtest_CXXFLAGS)
testAll_LDFLAGS =	$(gtest_LDFLAGS)
testAll_LDADD = $(gtest_LIBS) ../libkemin.la



