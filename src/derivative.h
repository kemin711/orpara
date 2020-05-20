#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <map>
#include <vector>
#include <exception>
#include <algorithm>

using namespace std;

/**
 * @return the derivative of 1st, 2nd order as a table of two number
 */
template<class T>
vector<tuple<T, double, double>> computeDerivative(const map<T, double> &data) {
   vector<pair<T,double>> input;
   copy(data.begin(), data.end(), back_inserter(input));
   return computeDerivative(input);
}
/**
 * @return 1st and 2nd order derivative
 */
template<class T>
vector<tuple<T, double, double>> computeDerivative(const vector<pair<T, double>> &input) {
   if (input.size() <= 3) {
      throw runtime_error(string(__FILE__) + ":ERROR: data size must be larger than 2");
   }
   // first and second derivative
   vector<tuple<T, double, double>> res(input.size());
   double d1, d2; 
   T dx, d2x;
   d1=(input[1].second - input[0].second)/(input[1].first - input[0].first);
   res[0] = make_tuple(input[0].first, d1, 0.0f);
   for (int i=1; i+1 < (int)input.size(); ++i) {
      d2x=input[i+1].first - input[i-1].first;
      d1 = (input[i+1].second - input[i-1].second)/d2x;
      d2 = 4*(input[i+1].second - 2*input[i].second + input[i-1].second)/(d2x*d2x);
      res[i] = make_tuple(input[i].first, d1, d2);
   }
   // compute derivative of the first data point
   get<2>(res[0]) = (get<1>(res[1]) - get<1>(res[0]))/(input[1].first - input[0].first);
   // last 1st derivative
   dx = input.back().first - input[input.size()-2].first;
   d1=(input.back().second-input[input.size()-2].second)/dx;
   d2=(d1-get<1>(res[res.size()-2]))/dx;
   res[res.size()-1] = make_tuple(input.back().first, d1, d2);
   return res;
}

#endif
