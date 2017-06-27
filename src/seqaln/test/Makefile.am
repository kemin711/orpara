TESTS = testall
check_PROGRAMS = testall 
testall_SOURCES = main_testall.cpp dynalnt_test.cpp 
testall_CXXFLAGS = -I.. -I../.. $(gtest_CXXFLAGS)
testall_LDFLAGS = $(gtest_LDFLAGS)
testall_LDADD = $(gtest_LIBS) ../../libkemin.la
