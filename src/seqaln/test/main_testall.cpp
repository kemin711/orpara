#include <iostream>
#include <gtest/gtest.h>
#include <bioseq.h>
#include <dynalnt.h>

using namespace std;
using namespace orpara;

int main(int argc, char* argv[]) {
   MatrixScoreMethod::setDefaultPath("../matrix");

   testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}

