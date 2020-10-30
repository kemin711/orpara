#include <iostream>
#include <gtest/gtest.h>
#include "quicksortorder.h"
#include <iterator>
#include <functional>
#include <algorithm>
#include <chrono>
#include <random>

using namespace orpara;
using namespace std;

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
   for (auto i=0; i<q_idx; ++i) {
      cout << input[i] << " ";
   }
   cout << "\nlarger portion:\n";
   for (int i=q_idx+1; i < (int)input.size(); ++i) {
      cout << input[i] << " ";
   }
   cout << endl;
   EXPECT_EQ(13, input[q_idx]);
   EXPECT_EQ(9, q_idx);
}
TEST(QuicksortorderTest, PartitionLargest) {
   vector<int> input={99, 1,2,1,2,3,8};
   int q_idx = partition<int>(input, 0, input.size()-1);
   cout << "pivot element: " << q_idx << " " << input[q_idx] << endl;
   cout << "smaller portion:\n";
   for (int i=0; i<q_idx; ++i) {
      cout << input[i] << " ";
   }
   cout << "\nlarger portion:\n";
   if (q_idx == int(input.size()-1)) {
      cout << "empty\n";
   }
   else {
      for (size_t i=q_idx+1; i<input.size(); ++i) {
         cout << input[i] << " ";
      }
   }
   cout << endl;
   EXPECT_EQ(99, input[q_idx]);
   EXPECT_EQ(6, q_idx);
}

TEST(QuicksortorderTest, PartitionSmallest) {
   vector<int> input={0, 1,2,1,1,2,3,8,9,9};
   int q_idx = partition<int>(input, 0, input.size()-1);
   cout << "pivot element: " << q_idx << " " << input[q_idx] << endl;
   cout << "smaller portion:\n";
   if (q_idx == 0) {
      cout << "empty\n";
   }
   else {
      for (int i=0; i<q_idx; ++i) {
         cout << input[i] << " ";
      }
   }
   cout << "\nlarger portion:\n";
   for (size_t i=q_idx+1; i<input.size(); ++i) {
      cout << input[i] << " ";
   }
   cout << endl;
   EXPECT_EQ(0, input[q_idx]);
   EXPECT_EQ(0, q_idx);
}

TEST(QuicksortorderTest, RandomPartitionVector) {
   using namespace std::placeholders;
   vector<int> input={13,19,9,35,12,8,7,4,11,2,6,21};
   cout << "Input vector\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   int q_idx = RandomPartition<int>()(input, 0, input.size()-1);
   cout << "Random pivot element: " << q_idx << " (" << input[q_idx] << ")\n";
   cout << "smaller portion:\n";
   for (int i=0; i<q_idx; ++i) {
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

// test the function version
TEST(QuicksortorderTest, RandomPartitionFunction) {
   using namespace std::placeholders;
   vector<int> input={13,19,9,35,12,8,7,4,11,2,6,21};
   cout << "Input vector\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   int q_idx = randomPartition<int>(input, 0, input.size()-1);
   cout << "Random pivot element: " << q_idx << " " << input[q_idx] << endl;
   cout << "smaller portion:\n";
   for (int i=0; i<q_idx; ++i) {
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
   cout << endl << endl;
   FindMedian<int> fm(input);
   // median input : {13,19,9,5,12,38,7,4,31,2,6,21};
   // after sort: {2,4,5,6,7,9,12,13,19,21,31,38};
   // 2 | 4 | 5 | 6 | 7 | 9 | 12 | 13 | 19 | 21 | 31 | 38 |
   int medianVal = fm.getMedian();
   cout << "Median value: " << medianVal << " should be 10 " << endl;
   ASSERT_EQ(medianVal, 10);
   cout << endl;
   fm(35);
   // after sort: {2,4,5,6,7,9,12,13,19,21,31,35,38};
   medianVal = fm.getMedian();
   cout << "Median value after adding 35: " << medianVal << " should be 12" << endl;
   ASSERT_EQ(medianVal, 12);
   cout << endl;
}

TEST(QuicksortorderTest, FindMedianNewdata) {
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
   cout << "=======================\n";
   cout << endl;
   vector<int> newinput={23, 5, 25, 0, 42, 49, 8, 10, 11, 7, 58, 27, 0, 40, 0};
   cout << "new input\n";
   copy(newinput.begin(), newinput.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   fm.setData(newinput);
   // median       0  1  2     4     6   M   8       10      12      14
   // after sort: {0, 0, 0, 5, 7, 8, 10, 11, 23, 25, 27, 40, 42, 49, 58};
   medianVal = fm.getMedian();
   cout << "Median value of new data: " << medianVal << " should be 11" << endl;
   ASSERT_EQ(medianVal, 11);
   cout << "=======================\n";
   newinput = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 40, 23, 58, 7, 11, 10, 8, 49, 42, 25, 27};
   cout << "new input with a lots of zeros\n";
   copy(newinput.begin(), newinput.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   fm.setData(newinput);
   medianVal = fm.getMedian();
   cout << "Median value of many zeros: " << medianVal << " should be 2" << endl;
   ASSERT_EQ(medianVal, 2);
   cout << "=======================\n";
   // has a lot of more zeros than the median
   newinput={23, 5, 25, 0, 42, 49, 8, 10, 11, 7, 58, 27, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   copy(newinput.begin(), newinput.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   fm.setData(newinput);
   medianVal = fm.getMedian();
   cout << "Median value of many zeros: " << medianVal << " should be 0" << endl;
   ASSERT_EQ(medianVal, 0);
   cout << "=======================\n\n";

}

// many largest value of the same
TEST(QuicksortorderTest, FindMedianTandemLarge) {
   vector<int> input={13,31, 31, 19,31, 31, 31, 9,31,31,5,12,7,4,31,2,6};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   vector<int> sorted=input;
   randomQuicksort<int>(sorted);
   cout << "sorted vector:\n";
   copy(sorted.begin(), sorted.end(), ostream_iterator<int>(cout, " | "));
   cout << endl << endl;
   
   FindMedian<int> fm(input);
   //0   2   4      7  8  9           12      14    16
   //2,4,5,6,7,9,12,13,19,31, 31, 31, 31, 31, 31,31,31
   int medianVal = fm.getMedian();
   cout << "Median value: " << medianVal << " should be 19 " << endl;
   ASSERT_EQ(medianVal, 19);

   input.push_back(31);
   cout << "new input\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   fm.setData(input);
   //2,4,5,6,7,9,12,13,19,31, 31, 31, 31, 31, 31,31,31, 31
   medianVal = fm.getMedian();
   cout << "Median value of new data: " << medianVal << " should be 25" << endl;
   ASSERT_EQ(medianVal, 25);

   input.push_back(31);
   cout << "new input\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   fm.setData(input);
   //2,4,5,6,7,9,12,13,19,31, 31, 31, 31, 31, 31,31,31, 31
   medianVal = fm.getMedian();
   cout << "Median value of new data: " << medianVal << " should be 31" << endl;
   ASSERT_EQ(medianVal, 31);
}
// many small a few large
TEST(QuicksortorderTest, FindMedianOneTwo) {
   vector<int> input={2,1,2,1,2,1,2,2,7,8,9};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   vector<int> sorted=input;
   randomQuicksort<int>(sorted);
   cout << "sorted vector:\n";
   copy(sorted.begin(), sorted.end(), ostream_iterator<int>(cout, " | "));
   cout << endl << endl;
   
   FindMedian<int> fm(input);
   // 1 2 3 4 5   7   9   11
   //{1,1,1,2,2,2,2,2,7,8,9};
   int medianVal = fm.getMedian();
   cout << "Median value (" << medianVal << ") should be 2 " << endl;
   ASSERT_EQ(medianVal, 2);
}
TEST(QuicksortorderTest, FindMedianBig) {
   vector<int> input={0,1,2,1,1,3,8,9,9,8,9,8,9,9,9,9,9,9,9,9};
   cout << "vector before sort\n";
   copy(input.begin(), input.end(), ostream_iterator<int>(cout, " | "));
   cout << endl;
   vector<int> sorted=input;
   randomQuicksort<int>(sorted);
   cout << "sorted vector:\n";
   copy(sorted.begin(), sorted.end(), ostream_iterator<int>(cout, " | "));
   cout << endl << endl;
   
   FindMedian<int> fm(input);
   // 1 2 3 4 5   7   9   11  13  15  17
   //
   //{0,1,1,1,2,3,8,8,8,9,9,9,9,9,9,9,9};
   int medianVal = fm.getMedian();
   cout << "Median value (" << medianVal << ") should be 9 " << endl;
   ASSERT_EQ(medianVal, 9);
}

TEST(QuicksortorderTest, FindMedianSpeed) {
   // test to see it is linear or n*log(n)
   random_device r;
   default_random_engine el(r());
   uniform_int_distribution<int> uniform_dist(1, 150);
   for (int n=100; n<10000000; n *= 2) {
      vector<int> data(n, 0);
      for (int i=0; i<n; ++i) {
         data[i] = uniform_dist(el);
      }
      chrono::time_point<chrono::steady_clock> t0 = chrono::steady_clock::now();
      FindMedian<int> fm(data);
      int medianV = fm.getMedian();
      chrono::steady_clock::duration tlen = chrono::steady_clock::now() - t0;
      cerr << "median " << medianV << " from " << n << " numbers took "
         << chrono::duration_cast<chrono::milliseconds>(tlen).count() << " milliseconds\n";
   }
}
      
