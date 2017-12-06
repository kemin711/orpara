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
   string rowstr="PR_RNA_NT1,PR_RNA_NT1,NA,NA,Pool24,NA,Pool24,NA,cfRNA_RT_20171201,PrediActV2,RNA,Predi_Run0067_NextSeq_171201_2,NA,5,1.38525";
   vector<string> row=splitQuoted(rowstr, '"', ',');
   cout << rowstr << endl;
   copy(row.begin(), row.end(), ostream_iterator<string>(cout, " | "));
}

TEST(StrformatTest, splitQuotedEmpty) {
   string rowstr=",Field2,PR_RNA_NT1,\"some thong to do with quote\",,,NA,Pool24,NA,Pool24,NA,cfRNA_RT_20171201,PrediActV2,RNA,Predi_Run0067_NextSeq_171201_2,NA,5,1.38525,LASTempty,";
   vector<string> row=splitQuoted(rowstr, '"', ',');
   cout << rowstr << endl;
   copy(row.begin(), row.end(), ostream_iterator<string>(cout, " | "));
}

int main(int argc, char* argv[]) {
   testing::InitGoogleTest(&argc, argv);
   // RUN_ALL_TEST macro from gtest
   return RUN_ALL_TESTS();
}

