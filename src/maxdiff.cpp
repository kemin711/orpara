#include <iostream>
#include "stddev.h"
#include <vector>

using namespace std;
using namespace orpara;

int main(int argc, char* argv[]) {
   vector<unsigned long> input{476232019, 487787336, 490498321, 500103870, 529907050, 540695814, 541056067, 546623159,
         552267810, 552818236, 557372854, 560951249, 565116326, 568780772, 574660216, 578322332,
         582110850, 595405819, 604857617, 618781588, 621258109, 623640992, 632250899, 632329067,
         636619716, 637425233, 639571346, 646165702, 650282086, 652572857, 652618365, 657486829,
         657637713, 664078742, 672483244, 673244078, 680649430, 692453644, 694445694, 694541912,
         708384169, 710398793, 720181320, 735814853, 15029402817, 15368636177, 27564355721, 28210994308};
   stddev avg1(input.front());
   stddev avg2(input.back());
   //  i      j
   //  ---=====
   int i=1;
   int j=input.size()-2;
   int plow, phigh;
   while (j > i) {
      double diffI1 = abs(avg1.getMean() - input[i]);
      double diffI2 = abs(avg2.getMean() - input[i]);
      if (diffI1 < diffI2) {
         avg1(input[i]);
         ++i;
         plow=i;
      }
      else { // all numbers from here on belong to the larger group
         // this is the end of the algorithm
         phigh=i;
         while (i < j) {
            avg2(input[i]);
            ++i;
         }
      }

      double diffJ1 = abs(avg1.getMean() - input[j]);
      double diffJ2 = abs(avg2.getMean() - input[j]);
      if (diffJ1 < diffJ2) {
         // j belong to the smaller group, all elements before j belong the the low group
         plow=j;
         while (j > i) {
            avg1(input[j]);
            --j;
         }
      }
      else {
         phigh = j;
         avg2(input[j]);
         --j;
      }
   }
   cout << " low at " << plow << " " << input[plow] << " hight at " 
      << phigh << " " << input[phigh] << " avg low high: "
      << avg1 << " " << avg2 << endl;
}

