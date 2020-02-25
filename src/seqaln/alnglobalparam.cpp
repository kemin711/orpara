#include "scorematrix.h"
#include "bioseq.h"
#include "dynalnt.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iterator>
//#include <stddev.h>

using namespace std;
using namespace orpara;

void scanGapParameter(const string &seq1, const string &seq2, int matchScore, int mismatchScore);
int multiplyRange(const pair<int,int> &r1, const pair<int,int> &r2);

/**
 * Probing the best parameters to use
 * match, mismatch, gapo gape (5, -4, -6, -4) indentity=0.782609
 * (4, -3, -5, -2) identity=0.782609 num nodes explored: 15532
 * for the test sequences.
 * best parameter: 4 3 0.782609 num nodes explored: 13799
 * (3, -3, -4, -3) 0.782609 num nodes explored: 13799
 */
int main(int argc, char* argv[]) {
   int i=1; 
   while (i < argc) {
      if (!strcmp(argv[i], "--help")) {
         cerr << "just test\n";
      }
      else {
         cerr << "your options: " << argv[i] << endl;
      }
      ++i;
   }
   string seq1="AGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTGAAAGGCAGTGGCTCAACCATTGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAGAAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAATACCGAT";
   string seq2="AGCGTTAATCGGAATTACTGGGCGTAAAGCGGGCGCAGACGGTTACTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCGTTTGAAACTGGGTGACTAGAGTATGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGAT";
   //string seq1="ACGGTTGC";
   //string seq2="AGCGTC";

   pair<int,int> matchRange={1, 12};
   pair<int,int> misRange={1,12};
   for (int i=matchRange.first; i<matchRange.second; ++i) {
      for (int j=misRange.first; j<misRange.second; ++j) {
         scanGapParameter(seq1, seq2, i, -j);
      }
   }

   return 0;
}

int multiplyRange(const pair<int,int> &r1, const pair<int,int> &r2) {
   return (r1.second-r1.first+1)*(r2.second-r2.first+1);
}

/** non random scan fast
 */
void scanGapParameter(const string &seq1, const string &seq2, int matchScore, int mismatchScore) {
   double maxidentity=0;
   //vector<int> numnodes;
   time_t timer1, timer2;
   time(&timer1);
   pair<int,int> gapOpenRange={1, 50};
   pair<int,int> gapExtRange={1,40};
   set<vector<int> > result; // go, ge, nnExplored
   for (int i=gapOpenRange.first; i<gapOpenRange.second; ++i) {
      int gapopen= -i;
      for (int j=gapExtRange.first; j<gapExtRange.second; ++j) {
         int gapext = -j;
         SimpleScoreMethod sm(matchScore, mismatchScore, -i, -j);
         DNA d1(seq1);
         DNA d2(seq2);
         Dynaln<SimpleScoreMethod> aligner(d1, d2);
         aligner.setMatrix(sm);
         aligner.runglobal();
         if (aligner.getIdentity()>maxidentity) {
            maxidentity = aligner.getIdentity();
            vector<int> tmp={gapopen, gapext};
            result.clear();
            result.insert(tmp);
         }
         else if (aligner.getIdentity() == maxidentity) {
            vector<int> tmp={gapopen, gapext};
            result.insert(tmp);
         }
      }
   }
   time(&timer2);
   double seconds = difftime(timer2, timer1);
   vector<string> header={"match", "mismatch", "gapopen", "gapext", 
      "identity", "time", "numTest"};
   copy(header.begin(), header.end(), ostream_iterator<string>(cout, "\t"));
   cout << endl;
   for (auto ptr=result.begin(); ptr != result.end(); ++ptr) {
      cout << matchScore << '\t' << mismatchScore << '\t';
      for (size_t j=0; j<(*ptr).size(); ++j) {
         cout << (*ptr)[j] << '\t';
      }
      cout << maxidentity << '\t' << seconds << '\t' 
         << multiplyRange(gapOpenRange, gapExtRange) << endl;
   }
}

