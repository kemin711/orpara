#include <iostream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>

/** testing the gsl library */
using namespace std;

int main(int argc, char* argv[]) {
   double x=5.0;
   double y=gsl_sf_bessel_J0(x);
   cout << y << endl;
   double zvalue[5]= {0.2, 0.75, -0.1, -0.3, 0.4};
   int nu=5;
   x=0;
   for (int i=0; i<nu; i++) {
      x += zvalue[i] * zvalue[i];
   }
   cout << "chi square tail " << x << " with " << nu << " degree of freedom "
      << gsl_cdf_chisq_P(x, nu) << endl;
      
   return 0;
}
