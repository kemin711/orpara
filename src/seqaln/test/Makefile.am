TESTS = testall testmatrix testkmerset testkmer
check_PROGRAMS = testall testmatrix testkmerset testkmer
testall_SOURCES = main_testall.cpp dynalnt_test.cpp 
testall_CXXFLAGS = -I.. -I../.. $(gtest_CXXFLAGS)
testall_LDFLAGS = $(gtest_LDFLAGS)
testall_LDADD = $(gtest_LIBS) ../../libkemin.la

testmatrix_SOURCES = scorematrix_test.cpp 
testmatrix_CXXFLAGS = -I.. -I../.. $(gtest_CXXFLAGS)
testmatrix_LDFLAGS = $(gtest_LDFLAGS)
testmatrix_LDADD = $(gtest_LIBS) ../../libkemin.la

testkmerset_SOURCES = kmerset_test.cpp
testkmerset_CXXFLAGS = -I.. -I../.. $(gtest_CXXFLAGS)
testkmerset_LDFLAGS = $(gtest_LDFLAGS)
testkmerset_LDADD = $(gtest_LIBS) ../../libkemin.la

testkmer_SOURCES = kmer_test.cpp
testkmer_CXXFLAGS = -I.. -I../.. $(gtest_CXXFLAGS)
testkmer_LDFLAGS = $(gtest_LDFLAGS)
testkmer_LDADD = $(gtest_LIBS) ../../libkemin.la

