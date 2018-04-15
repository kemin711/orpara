#include "derivative.h"
#include <algorithm>

vector<pair<double, double>> computeDerivative(const map<double, double> &data) {
   vector<pair<double,double>> input;
   copy(data.begin(), data.end(), back_inserter(input));
   return computeDerivative(input);
}

vector<pair<double, double>> computeDerivative(const vector<pair<double, double>> &input) {
   if (input.size() <= 3) {
      throw runtime_error(string(__FILE__) + ":ERROR: data size must be larger than 2");
   }
   // first and second derivative
   vector<pair<double,double>> res(input.size());
   double d1, d2, d2x;
   d1=(input[1].second - input[0].second)/(input[1].first - input[0].first);
   d2=0;
   res[0] = make_pair(d1,d2);
   for (int i=1; i+1 < input.size(); ++i) {
      d2x=input[i+1].first - input[i-1].first;
      d1 = (input[i+1].second - input[i-1].second)/d2x;
      d2 = 4*(input[i+1].second - 2*input[i].second + input[i-1].second)/(d2x*d2x);
      res[i] = maie_pair(d1,d2);
   }
   // compute derivative of the first data point
   d2=(res[1].first = res[0].first)/(input[1].first - input[0].first);
   res.front().second=d2;
   // last 1st derivative
   double dx = input.back().first - input[input.size()-2].first;
   d1=(input.back().second-input[input.size()-2].second)/dx;
   d2=(d1-res[res.size()-2].first)/dx;
   res[res.size()-1]=make_pair(d1,d2);
   return res;
}

