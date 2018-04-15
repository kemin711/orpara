#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <map>
#include <vector>
#include <exception>

using namespace std;

/**
 * @return the derivative of 1st, 2nd order as a table of two number
 */
vector<pair<double, double>> computeDerivative(const map<double, double> &data);
/**
 * @return 1st and 2nd order derivative
 */
vector<pair<double, double>> computeDerivative(const vector<pair<double, double>> &input);

#endif
