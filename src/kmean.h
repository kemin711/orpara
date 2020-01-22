#ifndef KMEAN_H 
#define KMEAN_H

#include <random>

using namespace std;

template<class T>
class Kmean {
   public:
      Kmean() : data(), cluster(), centriod(),
         rgen(std::random_device()), dice() {}
      Kmean(int nc, const vector<T>& point) 
         : data(point), cluster(nc), centriod(nc, array<double, T::dimension>),
           rgen(std::random_device()), dice(1, data.size()-2)
      {
         if (data.size() <= cluster.size()) {
            throw runtime_error("infufficient data");
         }
         set<int> used;
         used.insert(data.size()/2);
         for (int i=0; i< centriod.size(); ++i) {
            if (i == 0) {
               centriod[0]=&data[0];
            }
            else if (i == 1) {
               centriod[1] = &data[data.size()-1];
            }
            else if (i == 2) {
               centroid[2] = &data[data.size()/2];
            }
            else {
               int x = dice(rgen);
               while (used.find(x) != used.end()) 
                  x = dice(rgen);
               centriod[i] = data[x];
               used.insert(x);
            }
         }
      }
      ~Kmean() {
         //delete[] centriod;
      }
      void run();

   private:
      //int numcluster;
      vector<T> data;
      vector<vector<T*>> cluster;
      /**
       * dimention of data should be obtained
       * from T
       */
      vector<array<double, D>> centriod;
      std::mt19937 rgen;
      uniform_int_distribution<> dice;
      static const int maxiteration=50;
};

template<class T> void Kmean::run() {
   int minj;
   for (int i=0; i < data.size(); ++i) {
      double mindist=numeric_limits<double>::max();
      for (int j=0; j < centriod.size(); ++j) {
         double d = centroid[j]->distance(data[i]);
         if (d < mindist) {
            mindist=d; minj=j;
         }
      }
      cluster[minj].push_back(&data[i]);
   }
}


};

#endif
