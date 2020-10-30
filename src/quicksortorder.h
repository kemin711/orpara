#ifndef QUICKSORTORDER_H
#define QUICKSORTORDER_H

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <random>
#include <set>
#include <limits>

//#define DEBUGLOG

using namespace std;

namespace orpara {

/**
 * This function work on a range [p, r] p<r.
 *
 * @param arr is the input vector.
 * @param p is the starting index of the vector. 
 *       [0--L] where L is the last index of the vector.
 * @param r is the ending index of the vector.
 * @return the pivot index q such that
 *   all elements in [p, q-1] are no greater than arr[q]
 *   and all element in [q+1, r] are no smaller than arr[q]
 */
template<class T>
int partition(vector<T> &arr, int p, int r) {
   int x = arr[p];
   int i=p+1;
   int j=r;
#ifdef DEBUGLOG
   cout << __func__ << ":" << __LINE__ << ": Before partition: [" << p << ", " << r << "]\n";
   copy(arr.begin(), arr.end(), ostream_iterator<T>(cout, " | "));
   cout << endl;
   copy(arr.begin()+p, arr.begin()+r+1, ostream_iterator<T>(cout, " | "));
   cout << endl;
#endif
   while (true) {
#ifdef DEBUGLOG
      cout << " i=" << i << "(" << arr[i] << "), "
         << " j=" << j << "(" << arr[j] << ")\n";
#endif
      while (arr[j] > x && j > p) --j;
      while (arr[i] <= x && i < r) ++i;
#ifdef DEBUGLOG
      cout << "after walk i,j: " << i << ", " << j << endl;
#endif
      if (i < j) {
         swap(arr[i], arr[j]);
      }
      /* add this cause it to fail
      else if (i == j) {
         return j;
      }
      */
      else { // p could be the same as j
         swap(arr[p], arr[j]);
#ifdef DEBUGLOG
         cout << __func__ << ":" << __LINE__ << ": After partition: [" 
            << p << ", " << r << "]\n";
         copy(arr.begin(), arr.end(), ostream_iterator<T>(cout, " | "));
         cout << endl;
         cout << "sub range\n";
         copy(arr.begin()+p, arr.begin()+r+1, ostream_iterator<T>(cout, " | "));
         cout << endl;
         cout << "pivot idx: " << j << endl;
#endif
         return j;
      }
   }
}

template<class T>
class RandomPartition {
   public:
      static mt19937 generator_QSO;

      RandomPartition() { }
      int operator()(vector<T>& arr, int p, int r) {
         uniform_int_distribution<int> unif_dist(p, r);
         int i = unif_dist(generator_QSO);
         swap(arr[i], arr[p]);
         return partition<T>(arr, p, r);
      }
};

template<class T>
mt19937 RandomPartition<T>::generator_QSO=mt19937(random_device()());


template<class T>
int randomPartition(vector<T> &arr, int p, int r) {
   //default_random_engine generator(random_device()());
   std::random_device rd;
   std::mt19937 gen(rd());
   uniform_int_distribution<int> unif_dist(p, r);
   int i = unif_dist(gen); // the first number is always p
   //int i = unif_dist(generator);
   //int i = unif_dist(generator_QSO);
   //cout << "randome idx: " << i << " in [" << p << ", " << r << "]\n";
   swap(arr[i], arr[p]);
   //cout << "after swap\n";
   //copy(arr.begin(), arr.end(), ostream_iterator<T>(cout, " | "));
   //cout << endl;
   return partition<T>(arr, p, r);
}

template<class T>
void quicksort(vector<T> &arr, int p, int r) {
   if (p < r) {
      if (r-p == 1) {
         //cout << __func__ << ":solving base case\n";
         if (arr[p] <= arr[r]) return;
         else swap(arr[p], arr[r]);
      }
      else {
         int q = partition(arr, p, r);
         quicksort<T>(arr, p, q-1);
         quicksort<T>(arr, q+1, r);
      }
   }
}

template<class T>
void quicksort(vector<T> &arr) {
   quicksort(arr, 0, arr.size()-1);
}

template<class T>
void randomQuicksort(vector<T> &arr, int p, int r) {
   if (p < r) {
      if (r-p == 1) {
         //cout << __func__ << ":solving base case\n";
         if (arr[p] <= arr[r]) return;
         else swap(arr[p], arr[r]);
      }
      else {
         //int q = randomPartition(arr, p, r);
         int q = RandomPartition<T>()(arr, p, r);
         randomQuicksort<T>(arr, p, q-1);
         randomQuicksort<T>(arr, q+1, r);
      }
   }
}

template<class T>
void randomQuicksort(vector<T> &arr) {
   randomQuicksort(arr, 0, arr.size()-1);
}

/**
 * select the ith smallest number from array's
 * [p, r] range
 * @param arr input vector.
 * @param i is the ith smallest number.
 * @return the ith smallest value in the array.
 */
template<class T>
int randomSelect(vector<T> &arr, int p, int r, int i) {
   //cout << __func__ << " picking " << i << "th element from "
   //   << " [" << p << ", " << r << "]\n";
   if (p == r) return arr[p];
   int q = RandomPartition<T>()(arr, p, r);
   int k = q-p+1;
   if (i <= k) 
      return randomSelect<T>(arr, p, q, i);
   else 
      return randomSelect<T>(arr, q+1, r, i-k);
}

/**
 * None recursive version not limited by the 
 * calling stack size.
 */
template<class T>
int randomSelectIterative(vector<T> &arr, int p, int r, int i) {
   while (p != r) {
      int q = RandomPartition<T>()(arr, p, r);
      int k = q-p+1;
      if (i <= k) r = q;
      else {
         p = q+1;
         i -= k;
      }
   }
   return arr[p];
}

/**
 * Algorithm object to calculate the median for
 * a given data.
 */
template<class T>
class FindMedian {
   public:
      /**
       * Default constructor with empty data
       */
      FindMedian() : data(), medianVal(numeric_limits<T>::max()), numuniq(-1) 
      { }
      /**
       * Initialize from a single datum.
       */
      FindMedian(const T& d)
         : data(1,d), medianVal(numeric_limits<T>::max()), numuniq(-1) 
      { }
      /**
       * @param n number of the same value d
       */
      FindMedian(const T& d, int n)
         : data(n,d), medianVal(numeric_limits<T>::max()), numuniq(-1) 
      { }
      /**
       * Constructor from a vector as input.
       */
      FindMedian(const vector<T>& arr) 
         : data(arr), medianVal(numeric_limits<T>::max()), numuniq(-1) 
      {}
      FindMedian(const FindMedian& other)
         : data(other.data), medianVal(other.medianVal), numuniq(other.numuniq)
      { }
      FindMedian(FindMedian&& other)
         : data(std::move(other.data)), 
           medianVal(other.medianVal), numuniq(other.numuniq)
      { }
      /**
       * Accumulate indidual data.
       */
      void operator()(T val) {
         data.push_back(val);
      }
      void operator()(T val, int n) {
         data.insert(data.cend(), n, val);
      }

      FindMedian& operator=(const FindMedian& other) {
         if (this != &other) {
            data=other.data;
            medianVal=other.medianVal;
            numuniq=other.numuniq;
         }
         return *this;
      }

      FindMedian& operator=(FindMedian&& other) {
         if (this != &other) {
            data=std::move(other.data);
            medianVal=other.medianVal;
            numuniq=other.numuniq;
         }
         return *this;
      }

      /**
       * Computation will alter the stats of 
       * the data member.
       */
      T getMedian() {
         /*
         if (medianVal != numeric_limits<T>::max()) {
#ifdef DEBUGLOG
            cout << "median value already computed: " << medianVal << endl;
#endif
            return medianVal;
         }
         */
         if (data.empty()) {
#ifdef DEBUGLOG
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
               << ":WARN medianval=" << medianVal << " although data is empty, please check your logic\n";
#endif
            return medianVal;
         }
         else if (data.size() == 1) {
            medianVal=data.front();
         }
         else if (data.size() == 2) {
            medianVal = (data[0]+data[1])/2;
         }
         else {
            int Im = data.size()/2;
            medianVal = randomselIterative(data, 0, data.size()-1, Im+1);
            if (data.size() % 2 == 0) {
               T medianAfter = randomselIterative(data, 0, data.size()-1, Im);
               medianVal = (medianAfter+medianVal)/2;
            }
         }
#ifdef DEBUGLOG
         cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << endl;
         show(cout);
#endif
         return medianVal;
      }

      /**
       * Pure function to compute the median value.
       */
      T computeMedian() {
         if (data.empty()) {
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
               << ":WARN data is empty, please check your logic\n";
            medianVal=numeric_limits<T>::max();
         }
         else if (data.size() == 1) { // single value
            medianVal=data.front();
         }
         else if (data.size() == 2) {
            medianVal = (data[0]+data[1])/2;
         }
         else {
            int Im = data.size()/2;
            medianVal = randomselIterative(data, 0, data.size()-1, Im+1);
            if (data.size() % 2 == 0) {
               T medianAfter = randomselIterative(data, 0, data.size()-1, Im);
               medianVal = (medianAfter+medianVal)/2;
            }
         }
         return medianVal;
      }

      T getMedianAndClear() {
         T tmp = getMedian();
         data.clear();
         return tmp;
      }

      void setData(const vector<T>& newdata) { 
         data = newdata; 
         numuniq=-1;
         medianVal=numeric_limits<T>::max();
      }

      /**
       * Clear the intermediate data 
       * and set the final result to 
       * default impossible numbers.
       */
      void clear() { 
         data.clear(); 
         medianVal=numeric_limits<T>::max();
         numuniq=-1;
      }

      bool empty() const { return data.empty(); }
      int getN() const { return data.size(); }

      int getUniqueCount() const {
         if (data.empty()) {
            if (numuniq != -1) return numuniq;
#ifdef DEBUGLOG
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << ":WARN empty data, check logic\n";
#endif
            return -1;
         }
         //if (numuniq == -1) {
            set<T> tmp(data.begin(), data.end());
            numuniq=tmp.size();
#ifdef DEBUGLOG
            cerr << __FILE__ << ":" << __LINE__ << " data size: " 
               << data.size() << " unique: " << numuniq << endl;
#endif
         //}
         return numuniq;
      }

      /**
       * Compute the unique count and return the value.
       */
      int computeUniqueCount() const {
         if (data.empty()) {
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << ":WARN empty data, check logic\n";
            numuniq=-1;
         }
         else {
            set<T> tmp(data.begin(), data.end());
            numuniq=tmp.size();
         }
         return numuniq;
      }

      pair<T, int> getMedianUnique() const {
         return make_pair(getMedian(), getUniqueCount());
      }

      /**
       * This method can only be called once!
       *
       * Compute and store median and numberOfUnique in this
       *  object. Clear the storage, then return the computed values.
       * @return the [median, uniqueCount] pair
       *    then clear the data member to save storage.
       */
      pair<T, int> getMedianUniqueAndClear() {
         computeMedian();
         computeUniqueCount();
         data.clear();
         return make_pair(medianVal, numuniq);
      }

      pair<T, int> computeMedianUnique() {
         computeMedian();
         computeUniqueCount();
         return make_pair(medianVal, numuniq);
      }

      void showData(ostream &ous) const {
         if (data.empty()) {
            ous << "data empty\n";
         }
         else {
            for (const T& d : data) {
               ous << d << " | ";
            }
            ous << endl;
         }
      }

      void show(ostream& ous) const {
         ous << "FindMedian Object\n"
            << "data: ";
         copy(data.begin(), data.end(), ostream_iterator<T>(ous, ", "));
         //showData(ous);
         ous << "\nmedian=" << medianVal << " numUnique=" << numuniq
            << endl;
      }

      static void showArray(const vector<T>& arr, ostream& ous) {
         copy(arr.begin(), arr.end(), ostream_iterator<T>(ous, "|"));
      }

   private:
      /**
       * Data member.
       */
      vector<T> data;
      /**
       * Median value from data. After obtaining this value
       * the data can be cleared.
       */
      mutable T medianVal;
      /**
       * Number of unique values
       */
      mutable int numuniq;

      static int randompart(vector<T> &arr, int p, int r) {
         uniform_int_distribution<int> unif_dist(p, r);
         int i = unif_dist(rand_engine);
#ifdef DEBUGLOG
         cout << __func__ << ": random idx: " << i << " in [" << p << ", " << r << "]\n";
#endif
         swap(arr[i], arr[p]);
         return partition<T>(arr, p, r);
      }

      /**
       * @param i is the ith starting from 1st, 2nd, ...
       * @param arr input array
       * @param p start of the array
       * @param r end of the array [p, r] close range.
       *
       * This is a helper function
       */
      static int randomselIterative(vector<T> &arr, int p, int r, int i) {
#ifdef DEBUGLOG
         cout << __func__ << ": input selecting " << i << "th\n";
         showArray(arr, cout);
#endif
         while (p != r) {
#ifdef DEBUGLOG
            //cout << endl << __func__ << ": p=" << p << ", r=" << r << " i=" << i << endl;
#endif
            if (r-p == 1) {
               if (arr[r] < arr[p]) {
                  swap(arr[p], arr[r]);
               }
               if (i == 1) return arr[p];
               if (i == 2) return arr[r];
               cerr << "ERROR: i " << i << " can only be 1 or 2\n";
               exit(1);
            }
            int q = randompart(arr, p, r);
            if (q == p) {
               if (i == 1) return arr[p];
               ++p;
               --i;
               continue;
            }
            else if (q == r) {
               if (i == r-p+1)  {
                  return arr[r];
               }
               --r;
               continue;
            }
            int k = q-p+1;
            if ( k == i) {
               return arr[q];
            }
            if (i < k) {
               r = q-1;
            }
            else {
               p = q+1;
               i -= k;
            }
         }
#ifdef DEBUGLOG
         cout << __func__ << ":" << __LINE__
            << " selected value: " << arr[p] << endl;
#endif
         return arr[p];
      }

      static mt19937 rand_engine;
};

template<class T>
mt19937 FindMedian<T>::rand_engine = mt19937(random_device()());

} // end of orpara namespace

#endif
