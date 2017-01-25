#ifndef KMERHELPER_H
#define KMERHELPER_H

// (c) 2012 Kemin Zhou at orpara.com

#include <vector>
#include <string>
#include <iostream>
#include <set>

using namespace std;

namespace orpara {
/**
 *  convert sequence to an array indexed by the
 *  integer representation of the word
 *  @param s input string
 *  @param w word size
 */
vector<int> hashArray(const string &s, const int w);
/**
 * convert base to integer
 */
int b2i(char ch);
/**
 * helper function conversting k into mask for 
 * making the hash value of kmer.
 */

void check_pair(set<int> &loc);


unsigned int computeMask(int mersz);

/**
 * Used by kmercount
 */
double computeD2(const vector<double> &f1, const vector<double> &f2);
/**
 * remove numbers that don't have close neighbors
 */
void trimLonely(set<int> &loc);
/**
 * combine nearby numbers into large ranges.
 */
vector<pair<int,int> > combineRanges(set<int> &loc, int d);

vector<pair<int,int> > combineStemRanges(set<int> &loc);


pair<int,int> longestRegion(const vector<pair<int,int> > &reg);
pair<int,int> nearestRegion(const vector<pair<int,int> > &reg, int point);
}
#endif
