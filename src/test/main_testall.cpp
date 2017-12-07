#include <iostream>
#include <gtest/gtest.h>
#include <bioseq.h>
#include <string>
#include "../strformat.h"

using namespace std;
using namespace orpara;

TEST(GlobalFunctionTest, translate) {
   string peptideseq, nuclseq;
   nuclseq="ATGACGGCAAAATTTCCCTGA";
   translate(peptideseq, nuclseq, 1);
   ASSERT_EQ("MTAKFP*", peptideseq);
}

TEST(StrformatTest, splitQuotedNarrow) {
   string rowstr="XRNA_NT1,PR_NT1,NA,NA,Pool24,NA,Pool24,NA,RNA_RT_20171201,tacctV2,RNA,Run0067_Seq_171201_2,NA,5,1.38525";
   vector<string> row=splitQuoted(rowstr, '"', ',');
   cout << rowstr << endl;
   copy(row.begin(), row.end(), ostream_iterator<string>(cout, " | "));
}

TEST(StrformatTest, splitQuotedEmpty) {
   string rowstr=",Field2,PR_RNA_NT1,\"some thong to do with quote\",,,NA,Pool24,NA,Pool24,NA,RNA_RT_20171201,XYX,RNA,Run0067_xeq_171201_2,NA,5,1.38525,LASTempty,";
   vector<string> row=splitQuoted(rowstr, '"', ',');
   cout << rowstr << endl;
   copy(row.begin(), row.end(), ostream_iterator<string>(cout, " | "));
}

TEST(StrformatTest, splitQuotedTail) {
   string input="Run0066_exteq_171104_1,,";
   vector<string> row=splitQuoted(input, '"', ',');
   cout << input << endl;
   copy(row.begin(), row.end(), ostream_iterator<string>(cout, " | "));
   ASSERT_EQ(row.size(), 3);
}
   

int main(int argc, char* argv[]) {
   testing::InitGoogleTest(&argc, argv);
   // RUN_ALL_TEST macro from gtest
   return RUN_ALL_TESTS();
}

