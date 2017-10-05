#include <getest/gtest.h>
#include "kmerset.h"

using namespace std;
using namespace orpara;

class KmersetTest : public testing::Test {
   protected:
      KmersetTest() :
         i7adapter("AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"),
         i5adapter("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"),
         adp(), kmst()
      { }

      virtual void SetUp() {
         adp.setAdapterFile("/localdrive/lbwfdata/Adapter/indexpool.txt");
         kmst.eat(i7adapter);
         kmst.eat(i5adapter);
         vector<string> barcodes = adp.getBarocodeVector();
         for (auto& bc : barcodes) {
            kmst.eat(bc);
         }
         kmst.addRC();
      }
      virtual void TearDown() { }

      // raw data
      string i7adapter;
      string i5adapter;
      AdapterPool adp;
      KmerSet<5> kmst;
};

TEST_F(KmersetTest, common) {
   // right adapter at 136
   string rightAdapter="CCCAGAGGCCTTCATGGAAGGAATATTCACTTCTAAAACAGACACATGGTAAGTCAGCCATCATCCTCCAGGTATCCCTGCAGCCATAAGGTGGTGCTCCTGGGCCAAAGGACTCTATACTCTAAGCCGGGAGCCCAGATCGGAAGAGCAC";
   int nc = kmst.common(rightAdapter);
   assert(nc > 1);
}
