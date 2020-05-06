#include <gtest/gtest.h>
#include "kmerbase.h"
#include <string>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <strformat.h>
#include <cstdint>

using namespace std;
using namespace orpara;

class KmerTest : public testing::Test {
   protected:
      KmerTest() :
         rawseq("AATGATACGGCGACCACCGAGATCTACACTGACAGCTATAACACTCTTTCCCTACACGACGCTCTTCCGATCT")
      { }

      virtual void SetUp() {
      }

      virtual void TearDown() { }

      // raw data
      string rawseq;
};

TEST_F(KmerTest, common) {
   string rightAdapter="CCCAGAGGCCTTCATGGAAGGAATATTCACTTCTAAAACAGACACATGGTAAGTCAGCCATCATCCTCCAGGTATCCCTGCAGCCATAAGGTGGTGCTCCTGGGCCAAAGGACTCTATACTCTAAGCCGGGAGCCCAGATCGGAAGAGCAC";
   string leftAdapter("TGCTCTTCCTATCTAAATTTGACAAAAGTATTCACTGTTCCATAATGAAGTTAATGTCTCCACCACTGGATTTCTCAGGAATCACTGACATAGGAGAAGTTTCCCAATTTCTGACCGAGGGAATCATCATGAAAGATTTTAGTCATCCCAA");
   uint64_t msk = KmerBase<3>::getMask();
   uint64_t l_msk = KmerBase<3>::getLeftMask();
   uint64_t r_msk = KmerBase<3>::getRightMask();
   cout << std::hex << msk << " " << l_msk << " " << r_msk << endl;
   ASSERT_TRUE(msk == 63);
}
