#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <cmath>
/** 
 * @Evolution
 * This package contains functions for generating
 * phylogenetic trees
 */

/** 
 * Estimate the number if substitutions per site
 * given the observed fraction of identical
 * residues between two sequences. Using the
 * Grishin (1995) method.
 *
 * @param q is the probility of the site
 *    being the same. Is is usually measured
 *    as the fraction of identical redidues
 * @param err is the error required for the 
 *    result. default 0.0000001
 * @return the number of amino acid substitutions per site
 */
double GrishinDistance(double q, double err=0.0000001) {
    double x=1-q;
    double xx=-1;
    while (std::abs(xx-x) > err) {
      xx=x;
      x -= (1+2*x)*(x*log(1+2*x)-2*x*x*q)/(2*x-(1+2*x)*log(1+2*x));
    }
    return x;
} 
#endif

