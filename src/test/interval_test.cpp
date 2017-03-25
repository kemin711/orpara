#include <iostream>
#include <gtest/gtest.h>
#include "interval.h"

using namespace orpara;

class IntervalTest : public testing::Test {
   protected:
      IntervalTest() 
         : empty(), foo(1, 7), bar(2, 9) { }

      virtual ~IntervalTest() { }
      virtual void SetUp() {
      }

      virtual void TearDown() {
      }

      Interval empty;
      Interval foo;
      Interval bar;
};


TEST_F(IntervalTest, Overlap) {
   int olp = foo.overlap(bar);
   ASSERT_EQ(olp, 6);
}

TEST_F(IntervalTest, Extend) {
   cout << foo << " " << bar << endl;
   int olp = foo.extend(bar);
   cout << " overlap " << olp << endl;
   cout << "after extend: " << foo << endl;
   ASSERT_EQ(olp, 6);
}

TEST_F(IntervalTest, Empty) {
   ASSERT_TRUE(empty.isNull());
}
