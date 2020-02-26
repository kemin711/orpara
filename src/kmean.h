#ifndef KMEAN_H 
#define KMEAN_H

#include <random>
#include <vector>
#include <list>
#include <iostream>
#include "stddev.h"

using namespace std;

namespace orpara {
/**
 * T could be a tuple
 */
template<class T> class Kmean {
   public:
      Kmean() : data(), cluster(), centroid(),
         rdevice(), rgen(std::random_device()), dice() {}

      /**
       * Constructor by making a copy of point.
       * This may be expensive. My use move to
       * move data into this object.
       */
      Kmean(int nc, const vector<T>& point) 
         : data(point), cluster(nc), 
            centroid(nc, nullptr),
            //new_centriod(nc, new double[T::dimension]),
            rdevice(),
           rgen(rdevice()), dice(1, data.size()-2)
      {
         init();
      }

      Kmean(int nc,  vector<T>&& point) 
         : data(std::move(point)), cluster(nc), 
            centroid(nc, nullptr),
            //new_centriod(nc, new double[T::dimension]),
           rgen(std::random_device()), dice(1, data.size()-2)
      {
         init();
      }

      ~Kmean() {
         for (size_t i=0; i<centroid.size(); ++i) {
            delete[] centroid[i];
            //delete[] new_centriod[i];
         }
      }

      void init() {
         if (cluster.size() < 2) {
            throw runtime_error("must have 2 or more cluters");
         }
         if (data.size() <= cluster.size()+3) {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__)
                 +  ":ERROR infufficient data");
         }
         for (auto& c : centroid) {
            c= new double[T::dimension];
         }
         data[0].assign(centroid[0]);
         data[data.size()-1].assign(centroid[1]);
         set<int> used;
         for (size_t i=2; i< centroid.size(); ++i) {
            if (i == 2) {
               data[data.size()/2].assign(centroid[2]);
               used.insert(data.size()/2);
            }
            else {
               int x = dice(rgen);
               while (used.find(x) != used.end()) 
                  x = dice(rgen);
               data[x].assign(centroid[i]);
               used.insert(x);
            }
         }
         //cout << "after init centroids\n";
         //showCentroid(cout);
         //cout << endl;
      }
      /**
       * Put data points to their nearest centroid
       */
      void divide();
      /**
       * The centriold will be updated with the mean
       * from each cluster.
       * @return the square difference between the old
       *  and tne newly computed centriods.
       */
      double update();
      /**
       * execute the algorithm
       */
      void run() {
         int numtry=0;
         double diff=numeric_limits<double>::max();
         while (numtry < maxiteration && diff > 0.01) {
            divide();
            diff = update();
            //cout << "total cnetroid move: " << diff << endl;
            ++numtry;
         }
         //cout << "at end of iteration numtry=" << numtry << endl;
      }

      list<vector<T>> getCluster() const {
         list<vector<T>> res;
         for (size_t c=0; c<cluster.size(); ++c) {
            vector<T> tmp; tmp.reserve(cluster[c].size());
            for (T* p : cluster[c]) {
               tmp.push_back(*p);
            }
            res.push_back(std::move(tmp));
         }
         return res;
      }

      /**
       * Debug function to show the final results
       */
      void showCluster(ostream& ous) const;
      void showCentroid(ostream& ous) const;

   private:
      //int numcluster;
      vector<T> data;
      /**
       * Number of cluster identical to number of centroid
       */
      vector<vector<T*>> cluster;
      /**
       * dimention of data should be obtained
       * from T
       * could use 1-dimentinal array to similar
       * 2-D array to eliminate the vector element.
       */
      //vector<double*> centroid;
      vector<double*> centroid;
      //vector<double*> new_centriod;
      std::random_device rdevice;
      std::mt19937 rgen;
      uniform_int_distribution<> dice;
      static const int maxiteration=100;
};

template<class T> void Kmean<T>::divide() {
   int minj;
   // clear each cluster
   for (auto& c : cluster) {
      c.clear();
   }
   for (size_t i=0; i < data.size(); ++i) {
      double mindist=numeric_limits<double>::max();
      for (size_t j=0; j < centroid.size(); ++j) {
         double d = data[i].distance(centroid[j]);
         if (d < mindist) {
            mindist=d; minj=j;
         }
      }
      cluster[minj].push_back(&data[i]);
   }
}

template<class T> double Kmean<T>::update() {
   double point[T::dimension];
   double sumDiff=0;
   for (size_t c=0; c<centroid.size(); ++c) {
      Stddev<T::dimension> avgCenter;
      for (const T* ptr : cluster[c]) {
         ptr->assign(point);
         avgCenter(point);
      }
      double diff=0;
      for (int i=0; i<T::dimension; ++i) {
         double d= avgCenter.getMean(i) - centroid[c][i];
         diff += d*d;
         // update centroid
         centroid[c][i] = avgCenter.getMean(i);
      }
      //cout << "centroid " << c << " squaire of diff: " << diff << endl;
      sumDiff += diff;
   }
   return sumDiff;
}

template<class T> void Kmean<T>::showCluster(ostream& ous) const {
   for (int c=0; c<cluster.size(); ++c) {
      ous << string(70, 'C') << endl;
      ous << "cluster: " << c << endl;
      for (T* p : cluster[c]) {
         p->show(ous);
         ous << endl;
      }
      ous << string(70, '=') << endl;
   }
   ous << endl;
}

template<class T> void Kmean<T>::showCentroid(ostream& ous) const {
   for (int i=0; i < centroid.size(); ++i) {
      ous << "centroid " << i << ": ";
      for (int j=0; j<T::dimension; ++j) {
         ous << centroid[i][j] << ", ";
      }
      ous << endl;
   }
}

} // end of name space orpara

#endif
