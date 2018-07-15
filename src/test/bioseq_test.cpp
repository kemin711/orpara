#include <iostream>
#include <gtest/gtest.h>
#include "bioseq.h"

using namespace orpara;

/**
 * By the convention proposed by the original author
 * Fixture should be named after OriginalClassNameTest,
 * such as FooTest
 */
class BioseqTest : public testing::Test {
   protected:
      BioseqTest() 
         : bsq_empty(), bsq_filled("seqname", "ACGTACGTAACTCAACTCCCCCGGGGGTACAACG") { }

      virtual void SetUp() {
      }

      virtual void TearDown() {
      }

      bioseq bsq_empty;
      bioseq bsq_filled;
};

TEST_F(BioseqTest, getname) {
   const string seqname = bsq_filled.getName();
   ASSERT_EQ("seqname", seqname);
}

TEST_F(BioseqTest, subseq) {
   bioseq sub = bsq_filled.subseq(3, 12);
   ASSERT_EQ(sub, bioseq("GTACGTAACT"));
}

TEST_F(BioseqTest, empty) {
   ASSERT_EQ(true, bsq_empty.empty());
}
