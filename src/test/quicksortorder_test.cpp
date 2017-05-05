#include <iostream>
#include <gtest/gtest.h>
#include "quicksortorder.h"
#include <iterator>
#include <functional>
#include <algorithm>

using namespace orpara;

TEST(QuicksortorderTest, RandomNumberGenerator) {
   default_random_engine generator;
   std::random_device rd;
   std::mt19937 gen(rd());
   cout << "random number from 0 to 100 ...\n";
   uniform_int_distribution<int> unif_dist(0, 100);
   //int i = unif_dist(generator);
   int i = unif_dist(gen);
   cout << "first value: " << i << endl;
   for (size_t j=0; j<20; ++j) {
      //cout << "round: " << j << " " << unif_dist(generator) << endl;
      cout << "round: " << j << " " << unif_dist(gen) << endl;
   }
   //i = unif_dist(generator);
   i = unif_dist(gen);
   cout << "last i value: " << i << endl;

   ASSERT_TRUE(i>=0 && i<=100);
}

TEST(QuicksortorderTest, PartitionVector) {
   vector<int> input={13,19,9,5,12,8,7,4,11,2,6,21};
   int q_idx = partition<int>(input, 0, input.size()-1);
   cout << "pivot element: " << q_idx << " " << input[q_idx] << endl;
   cout << "smaller portion:\n";
   for (size_t i=0; i<q_idx; ++i) {
      cout << input[i] << " ";
   }
   cout << "\nlarger portion:\n";
   for (size_t i=q_idx+1; i<input.size(); ++i) {
      cout << input[i] << " ";
   }
   cout << endl;
   EXPECT_EQ(13, input[q_idx]);
   EXPECT_EQ(9, q_idx);
}

TEST(QuicksortorderTest, RandomPartitionVector) {
   using namespace std::placeholders;
   vector<int> input={13,19,9,35,12,8,7,4,11,2,6,21};
   cout << "Input vector\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   int q_idx = randomPartition<int>(input, 0, input.size()-1);
   cout << "Random pivot element: " << q_idx << " " << input[q_idx] << endl;
   cout << "smaller portion:\n";
   for (size_t i=0; i<q_idx; ++i) {
      cout << input[i] << " ";
   }
   cout << "\nlarger portion:\n";
   for (size_t i=q_idx+1; i<input.size(); ++i) {
      cout << input[i] << " ";
   }
   cout << endl;
   int pivotVal = input[q_idx];
   ASSERT_TRUE(
         any_of(input.begin(), input.end(), bind(equal_to<int>(), _1, pivotVal) )
         );
}

TEST(QuicksortorderTest, QuicksortVector) {
   vector<int> input={13,19,9,5,12,38,7,4,11,2,6,21};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   quicksort<int>(input);
   cout << "vector after sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   vector<int> sorted={2,4,5,6,7,9,11,12,13,19,21,38};
   ASSERT_EQ(input, sorted);
}

TEST(QuicksortorderTest, RandomQuicksortVector) {
   vector<int> input={13,19,9,5,12,38,7,4,11,2,6,21};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   randomQuicksort<int>(input);
   cout << "vector after sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   vector<int> sorted={2,4,5,6,7,9,11,12,13,19,21,38};
   ASSERT_EQ(input, sorted);
}

TEST(QuicksortorderTest, RandomSelectVector) {
   vector<int> input={13,19,9,5,12,38,7,4,31,2,6,21};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   // 4th smallest: {13,19,9,5,12,38,7,4,31,2,6,21};
   // {2,4,5,6,7,9,12,13,19,21,31,38};
   int ith4 = randomSelect<int>(input, 0, input.size()-1, 4);
   cout << "4th smallest element: " << ith4 << endl;
   ASSERT_EQ(ith4, 6);
   int ith6 = randomSelect<int>(input, 0, input.size()-1, 6);
   cout << "6th smallest element: " << ith6 << endl;
   ASSERT_EQ(ith6, 9);
}

TEST(QuicksortorderTest, RandomSelectIterativeVector) {
   vector<int> input={13,19,9,5,12,38,7,4,31,2,6,21};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   // 4th smallest: {13,19,9,5,12,38,7,4,31,2,6,21};
   // {2,4,5,6,7,9,12,13,19,21,31,38};
   int ith4 = randomSelectIterative<int>(input, 0, input.size()-1, 4);
   cout << "randomSelectIterative() 4th smallest element: " << ith4 << endl;
   ASSERT_EQ(ith4, 6);
}

TEST(QuicksortorderTest, FindMedianVector) {
   vector<int> input={13,19,9,5,12,38,7,4,31,2,6,21};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   vector<int> sorted=input;
   randomQuicksort<int>(sorted);
   cout << "sorted vector:\n";
   copy(sorted.begin(), sorted.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   FindMedian<int> fm(input);
   // median input : {13,19,9,5,12,38,7,4,31,2,6,21};
   // after sort: {2,4,5,6,7,9,12,13,19,21,31,38};
   // 2 | 4 | 5 | 6 | 7 | 9 | 12 | 13 | 19 | 21 | 31 | 38 |
   int medianVal = fm.getMedian();
   cout << "Median value: " << medianVal << " should be 10 " << endl;
   ASSERT_EQ(medianVal, 10);
   fm(35);
   // after sort: {2,4,5,6,7,9,12,13,19,21,31,35,38};
   medianVal = fm.getMedian();
   cout << "Median value after adding 35: " << medianVal << " should be 12" << endl;
   ASSERT_EQ(medianVal, 12);
}

