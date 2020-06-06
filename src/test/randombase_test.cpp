#include <iostream>
#include <gtest/gtest.h>
#include "randombase.h"

using namespace orpara;

class RandomBaseTest : public testing::Test {
   protected:
      RandomBaseTest() : rbengine(RandomBase::getInstance()) { }

      virtual ~RandomBaseTest() { }
      virtual void SetUp() {
      }

      virtual void TearDown() {
      }
      RandomBase& rbengine; 
};


TEST_F(RandomBaseTest, basen) {
   map<char,int> bcount;
   cout << "Random ACGT\n";
   for (auto i=0; i < 250; ++i) {
      cout << rbengine() << " ";
      ++bcount[rbengine()];
   }
   cout << endl;
   for (auto& p : bcount) {
      cout << p.first << " " << p.second << endl;
   }
   cout << endl;
   ASSERT_TRUE(abs(bcount['A'] - bcount['T']) < 70);
}

TEST_F(RandomBaseTest, basey) { // C or T
   map<char,int> bcount;
   cout << "Random CT\n";
   for (auto i=0; i < 850; ++i) {
      char x = rbengine('Y');
      cout << x << " ";
      ++bcount[x];
   }
   cout << endl;
   for (auto& p : bcount) {
      cout << p.first << " " << p.second << endl;
   }
   cout << endl;
   ASSERT_TRUE(abs(bcount['C'] - bcount['T']) < 78);
   ASSERT_TRUE(bcount.size() == 2);
}

TEST_F(RandomBaseTest, basew) { // C or T
   map<char,int> bcount;
   cout << "Random AT\n";
   for (auto i=0; i < 990; ++i) {
      char x = rbengine('W');
      cout << x << " ";
      ++bcount[x];
   }
   cout << endl;
   for (auto& p : bcount) {
      cout << p.first << " " << p.second << endl;
   }
   cout << endl;
   ASSERT_TRUE(abs(bcount['A'] - bcount['T']) < 78);
   ASSERT_TRUE(bcount.size() == 2);
}

TEST_F(RandomBaseTest, basev) { // A, C or G
   map<char,int> bcount;
   cout << "Random ACG\n";
   for (auto i=0; i < 950; ++i) {
      char x = rbengine('V');
      cout << x << " ";
      ++bcount[x];
   }
   cout << endl;
   for (auto& p : bcount) {
      cout << p.first << " " << p.second << endl;
   }
   cout << endl;
   ASSERT_TRUE(abs(bcount['C'] - bcount['G']) < 79);
   ASSERT_TRUE(bcount.size() == 3);
}

