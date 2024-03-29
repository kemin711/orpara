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

TEST_F(FastqTest, trimLowq) {
   // for low quality sequence, we need to use lower cutoff
   // 13 might be optimal
   Fastq fq("lowqualtity", "GATCATCCATACTTTTCTTTAGACGTTCTTCAGGATCAAGTTTTAATGATGCTGTGTGGCTGGATTTAAATTATCTTGAAGTTGCCAAGGGAGCTCAGTCTTGTGCTGCTCGCTTTACAGCTTTCCTCTATGCAGAACTCTATGCAGATAA", "AAAA6/A///////<//////E//////////EE//////<//<//AE/<E/E/A//EE/AEE///<////A/////E//E//E////EE//E////E///<E<E/6E6//////</</E/</////<//EE//E///<////E6/E/<//");
   cout << "original fastq:\n" << fq << endl;
   fq.trimLowq(7, 14);
   cout << "after trimLowq() fastq:\n" << fq << endl;
   ASSERT_GT(fq.length(), 15);
}

TEST_F(FastqTest, Gtail2) {
   Fastq fq("longGtail", "GCCTTTCCACCACCCCACCCTCGCCCTACCACCACCTTGGGGGCCGGGGGCGCGGGGGCGGGGGCGGGGGGGGGGGGGCGGGGGGGGCGGGGGGCCGCGGGGGGGGGGGGGGGGGGGGCGGGGGGGGCGGGGGGGGGGGGGGGCGG", "FFFFF6FF/FF/=FFF6FFA/FF=F/F/AF6F//F////F/////F/FAA/A//F/////FFAA/A/AAFA//AAF///FFF////F/AA/A/F//F//F/F///FAFA////F/F///A/F/F///////F///A//////////");
   cerr << "Fastq to be trimmed:\n" << fq << endl;
   cerr << "average quality=" << fq.getAverageQuality() << endl;
   if (fq.trimG()) {
      cerr << "G trimmed to\n" << fq << endl;
   }
}

TEST_F(FastqTest, trimN) {
   Fastq fq("withNtail", "ACATTTCAGACAGGAATTTTGTTCATTTTAATGAACTCCCACCATTCCAGCAGCTTTTTGTGATGATCCACATGTAATTGATTGTTCAAGAGATGCCTTATTTAACAAAGTACAGTGTACAAGCATACATAAGATTATGATNGNNNN", "AEEEAEEEE/EEEEEEEEEEEEEEEEEAEEEEEAAEEEEEEEEE/EEEEEEEEEAEAEEEEEEEEEEEEEEEEEEEAAEEEEEEAE<E<EEEEEAEEAEEEEAEEEEEEEEEEEE/EEE<EEE</AEE<EEEEAAAEEAEE#EEE#E");
   cerr << "Fastq to be trimmed of N:\n" << fq << endl;
   cerr << "average quality=" << fq.getAverageQuality() << endl;
   bool trimmed=false;
   if (fq.trimN()) {
      cerr << "N trimmed to\n" << fq << endl;
      trimmed = true;
   }
   ASSERT_TRUE(trimmed);
}
