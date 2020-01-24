#include <iostream>
#include <gtest/gtest.h>
#include <bioseq.h>
#include <string>
#include <fstream>
#include <list>

#include "../strformat.h"
#include "../derivative.h"
#include "../insertsortlist.h"


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

TEST(DerivativeTest, compute) {
   string inputFile="derivative_test.tab";
   vector<pair<int,double>> input;
   // first line is header x,y
   ifstream inf(inputFile);
   if (inf.fail()) {
      cerr << "Failed to open input file: " << inputFile << endl;
      exit(1);
   }
   string line;
   getline(inf, line); // remove header
   getline(inf, line);
   while (!inf.eof()) {
      string::size_type i = line.find('\t');
      input.push_back(make_pair(std::stoi(line.substr(0, i)), std::stof(line.substr(i+1))));
      getline(inf, line);
   }

   vector<tuple<int,double, double>> result=computeDerivative(input);
   //for (size_t i=0; i<input.size(); ++i) {
   //   cout << input[i].first << " " << input[i].second << " "
    //     << get<1>(result[i]) << " " << get<2>(result[i]) << endl;
   //}
}
 
TEST(InsertSortList, integerlist) {
   std::list<int> input1{3, 7, 12, 9, 4, 2, 1, 5, 7, 8};
   cout << "before sorting\n";
   for (int x : input1) {
      cout << x << ", ";
   }
   cout << endl;
   insertSortList(input1);
   for (int x : input1) {
      cout << x << ", ";
   }
   cout << endl;
}

TEST(InsertSortList, stringlist) {
   std::list<string> input1{"abcd", "xyz", "foo", "bar", "enjoy", "every", "1 moment", "live", "more", "888"};
   cout << "before sorting\n";
   for (auto& x : input1) {
      cout << x << ", ";
   }
   cout << endl;
   insertSortList(input1);
   for (auto& x : input1) {
      cout << x << ", ";
   }
   cout << endl;
   cout << "reverse sort\n";
   insertSortList(input1, [](const string& a, const string& b)->bool{ return a > b; });
   for (auto& x : input1) {
      cout << x << ", ";
   }
   cout << endl;
}

int main(int argc, char* argv[]) {
   testing::InitGoogleTest(&argc, argv);
   // RUN_ALL_TEST macro from gtest
   return RUN_ALL_TESTS();
}

