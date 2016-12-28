#include <iostream>
#include <gtest/gtest.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;

TEST(GlobalFunctionTest, translate) {
   string peptideseq, nuclseq;
   nuclseq="ATGACGGCAAAATTTCCCTGA";
   translate(peptideseq, nuclseq, 1);
   ASSERT_EQ("MTAKFP*", peptideseq);
}

int main(int argc, char* argv[]) {
   testing::InitGoogleTest(&argc, argv);
   // RUN_ALL_TEST macro from gtest
   return RUN_ALL_TESTS();
}

