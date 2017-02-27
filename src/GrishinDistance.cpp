#include <iostream>
#include <cmath>
#include "volution.h"

using namespace std;

//double numericFunc(double q, double error) {
//   double x=1-q;
//   double xx;
//   while (abs(xx-x) > error) {
//      xx=x;
//      x -= (1+2*x)*(x*log(1+2*x)-2*x*x*q)/(2*x-(1+2*x)*log(1+2*x));
//   }
//   return x;
//}

int main(int argc, char* argv[]) {
   double q;
   double err=0.0000001;

   if (argc>1) {
      q=atof(argv[1]);
      cerr << GrishinDistance(q, err) << endl;
   }
   else {
      for (double i=0.1; i<0.9999; i += 0.05) {
         cerr << i << " " << numericFunc(i, err) << endl;
      }
   }
   return 0;
}



