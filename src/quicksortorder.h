#ifndef QUICKSORTORDER_H
#define QUICKSORTORDER_H

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <random>
#include <set>

using namespace std;

namespace orpara{

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
   //cout << "Before partition: [" << p << ", " << r << "]\n";
   //copy(arr.begin(), arr.end(), ostream_iterator<T>(cout, " | "));
   //cout << endl;
   while (true) {
      //cout << " i,j " << i << "=" << arr[i] << ", " << j << "=" << arr[j] << endl;
      while (arr[j] > x && j > p) --j;
      while (arr[i] <= x && i < r) ++i;
      //cout << "after walk i,j: " << i << ", " << j << endl;
      if (i < j) {
         swap(arr[i], arr[j]);
      }
      /* add this cause it to fail
      else if (i == j) {
         return j;
      }
      */
      else {
         swap(arr[p], arr[j]);
         //cout << "After partition: [" << p << ", " << r << "]\n";
         //copy(arr.begin(), arr.end(), ostream_iterator<T>(cout, " | "));
         //cout << endl;
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

template<class T>
class FindMedian {
   public:
      FindMedian() : data() { }
      FindMedian(const vector<T>& arr) 
            : data(arr) {}
      void operator()(T val) {
         data.push_back(val);
      }

      T getMedian() {
         if (data.size() == 1) return data[0];
         if (data.size() == 2) return (data[0]+data[1])/2;
         int Im = data.size()/2;
         T medianVal = randomselIterative(data, 0, data.size()-1, Im+1);
         if (data.size() % 2 == 0) {
            T medianAfter = randomselIterative(data, 0, data.size()-1, Im);
            return (medianAfter + medianVal)/2;
         }
         return medianVal;
      }

      T getMedianAndClear() {
         T tmp = getMedian();
         clear();
         return tmp;
      }

      void setData(const vector<T>& newdata) { data = newdata; }

      void clear() { data.clear(); }
      bool empty() const { return data.empty(); }
      int getN() const { return data.size(); }
      int getUniqueCount() const {
         set<T> tmp;
         for (const T& d : data) {
            tmp.insert(d);
         }
         return tmp.size();
      }

      pair<T, int> getMedianUnique() {
         return make_pair(getMedian(), getUniqueCount());
      }

      pair<T, int> getMedianUniqueAndClear() {
         pair<T, int> tmp = make_pair(getMedian(), getUniqueCount());
         clear();
         return tmp;
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
         ous << "FindMedian: ";
         showData(ous);
      }

   private:
      vector<T> data;

      int randompart(vector<T> &arr, int p, int r) {
         uniform_int_distribution<int> unif_dist(p, r);
         int i = unif_dist(rand_engine);
         //cout << __func__ << ": randome idx: " << i << " in [" << p << ", " << r << "]\n";
         swap(arr[i], arr[p]);
         return partition<T>(arr, p, r);
      }

      int randomselIterative(vector<T> &arr, int p, int r, int i) {
         //showData(cout);
         while (p != r) {
            //cout << __func__ << ": p=" << p << ", r=" << r << " i=" << i << endl;
            int q = randompart(arr, p, r);
            //cout << "q returned by randompart() " << q << endl;
            int k = q-p+1;
            //cout << "k=" << k << endl;
            if (i <= k) {
               //cout << "i<=k" << endl;
               if (arr[p] == arr[r]) {
                  //cerr << "special case, median should be " << arr[r] << endl;
                  return arr[i-1];
                  //exit(1);
               }
               r = q;
            }
            else {
               //cout << "i>k" << endl;
               p = q+1;
               i -= k;
            }
         }
         //cout << "median value: " << arr[p] << endl;
         return arr[p];
      }
      static mt19937 rand_engine;
};

template<class T>
mt19937 FindMedian<T>::rand_engine = mt19937(random_device()());

} // end of orpara namespace

#endif