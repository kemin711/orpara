#include <gtest/gtest.h>
#include "kmerset.h"
#include <string>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <strformat.h>

using namespace std;
using namespace orpara;

class KmersetTest : public testing::Test {
   protected:
      KmersetTest() :
         i7adapter("AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"),
         i5adapter("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),
         adapterPoolFile("/localdrive/lbwfdata/Adapter/indexpool.txt"),
         kmst()
      { }

      virtual void SetUp() {
         kmst.eat(i7adapter);
         kmst.eat(i5adapter);
         set<string> barcodes = readBarcode();
         for (auto& bc : barcodes) {
            kmst.eat(bc);
         }
         kmst.addRC();
      }
      virtual void TearDown() { }

      set<string> readBarcode();

      // raw data
      string i7adapter;
      string i5adapter;
      //AdapterPool adp;
      string adapterPoolFile;
      //KmerSet<5> kmst;
      KmerSet<6> kmst;
};

// this test passed, got 102 common kmers
set<string> KmersetTest::readBarcode() {
   ifstream inf(adapterPoolFile);
   if (inf.fail()) {
      cerr << "Failed to read from " << adapterPoolFile << endl;
      exit(1);
   }
   set<string> res;
   string line;
   getline(inf, line);
   while (!inf.eof()) {
      vector<string> row=split(line, '\t');
      res.insert(row.begin()+1, row.end());
      getline(inf, line);
   }
   cout << "last line: " << line << endl;
   return res;
}


TEST_F(KmersetTest, common) {
   // right adapter at 136
   string rightAdapter="CCCAGAGGCCTTCATGGAAGGAATATTCACTTCTAAAACAGACACATGGTAAGTCAGCCATCATCCTCCAGGTATCCCTGCAGCCATAAGGTGGTGCTCCTGGGCCAAAGGACTCTATACTCTAAGCCGGGAGCCCAGATCGGAAGAGCAC";
   int nc = kmst.common(rightAdapter);
   cout << nc << " common kmers right\n";
   ASSERT_TRUE(nc > 1);
   string leftAdapter("TGCTCTTCCTATCTAAATTTGACAAAAGTATTCACTGTTCCATAATGAAGTTAATGTCTCCACCACTGGATTTCTCAGGAATCACTGACATAGGAGAAGTTTCCCAATTTCTGACCGAGGGAATCATCATGAAAGATTTTAGTCATCCCAA");
   nc = kmst.common(leftAdapter);
   cout << nc << " common kmers left\n";
   ASSERT_TRUE(nc > 1);
}
