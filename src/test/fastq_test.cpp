#include "fastq.h"
#include <string>
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>

using namespace std;
using namespace orpara;

class FastqTest : public testing::Test {
   protected:
      FastqTest() : emptyFq(), 
         grepeatFq("readWithManyG", "1:N:0:CTGATCGT+GCGCATAT", "TCTCTTCAGGGTCCCATGCTGGGCAGGAGGGGCCTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA////6FA/////F///6FFAFF/FF/////FFFFA/FAFFFAAFFFFA/FFFF/FFFF/FFFF/FF//F/FF///F//////AAFFFFF/////A/AAAFA/FA=//AA/AFA/FF")
          { }

      Fastq emptyFq;
      Fastq grepeatFq;
};

TEST_F(FastqTest, setter) {
   emptyFq.setName("foo");
   emptyFq.setTitle("Birth of a new fastq");
   emptyFq.setSequenceQuality("CTGTGGTGAGGGCTGAGGTGACCCTTGTCTCTGTGTTCTTGTCCCCCCCAGCTTGTGGAGCCTCTTACACCCAGTGGAGAAGCTCCCAACCAAGCTCTCTTGAGGATCTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCAG", "FFFFFFFFFFFFFFFFFFFFAFFFFFF/FFFFAFFFFFFFFFFFFFFFFFFFFF/FFFAFFFFFFFAF=FFFFFFFFFFAFFFFFFFFAFFFFFFFFFAAFFFFFFFFFAFFFFF/FFFFFAFFFFFFFFFFAFFFFFAFFFFFFFFAFFF");
   cout << "I made a fastq:\n" << emptyFq << endl;
   ASSERT_GT(emptyFq.length(), 0);
}

TEST_F(FastqTest, trimG) {
   cout << "before trim\n" << grepeatFq << endl;
   int oldlen=grepeatFq.length();
   bool trimmed = grepeatFq.trimG();
   if (trimmed) {
      cout << "after trim:\n" << grepeatFq << endl;
   }
   ASSERT_GT(oldlen, grepeatFq.length());
   oldlen=grepeatFq.length();
   trimmed = grepeatFq.trimLowq(5, 16);
   if (trimmed) {
      cout << "removed low qulaity region\n" << grepeatFq << endl;
      ASSERT_GT(oldlen, grepeatFq.length());
   }
}

