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
   cout << "mask: " << std::hex << msk << " " << l_msk << " " << r_msk << endl;

   uint64_t msk31 = KmerBase<31>::getMask();
   uint64_t l_msk31 = KmerBase<31>::getLeftMask();
   uint64_t r_msk31 = KmerBase<31>::getRightMask();
   cout << "31 bit mask: " << std::hex << msk31 << " " << l_msk31 << " " << r_msk31 << endl;

   uint64_t* kint = KmerBase<31>::kmerize(rawseq);
   size_t numKmer = rawseq.size()-30;
   for (size_t i=0; i<numKmer; ++i) {
      cout << std::dec << kint[i] << " ";
   }
   cout << endl;
   delete[] kint;

   ASSERT_TRUE(msk == 63);
}
