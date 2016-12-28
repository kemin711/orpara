#include "braboualn.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <stddev.h>

using namespace std;
void scanGapParameter(const string &seq1, const string &seq2, int matchScore, int mismatchScore);

/**
 * Probing the best parameters to use
 * match, mismatch, gapo gape (5, -4, -6, -4) indentity=0.782609
 * (4, -3, -5, -2) identity=0.782609 num nodes explored: 15532
 * for the test sequences.
 * best parameter: 4 3 0.782609 num nodes explored: 13799
 * (3, -3, -4, -3) 0.782609 num nodes explored: 13799
 */
int main(int argc, char* argv[]) {
   string seq1="AGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTGAAAGGCAGTGGCTCAACCATTGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAGAAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAATACCGAT";
   string seq2="AGCGTTAATCGGAATTACTGGGCGTAAAGCGGGCGCAGACGGTTACTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCGTTTGAAACTGGGTGACTAGAGTATGTCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGAT";
   //string seq1="ACGGTTGC";
   //string seq2="AGCGTC";

   vector<int> match={1,2,3,4,5,6,7,8,9,10};
   vector<int> mismatch={1,2,3,4,5,6,7,8,9,10};
   for (int i=0; i<match.size(); ++i) {
      for (int j=0; j<mismatch.size(); ++j) {
         scanGapParameter(seq1, seq2, match[i], -mismatch[j]);
      }
   }

   return 0;
}

void scanGapParameter(const string &seq1, const string &seq2, int matchScore, int mismatchScore) {
   srand(time(NULL));
   double maxidentity=0;
   vector<int> numnodes;
   time_t timer1, timer2;
   time(&timer1);
   int trial1=50;
   int trial2=50;
   set<vector<int> > result; // go, ge, nnExplored
   for (int i=0; i<trial1; ++i) {
      int gapopen= -(rand() % 22 + 1);
      for (int j=0; j<trial2; ++j) {
         int gapext = -(rand() % 17 + 1);
         Braboualn aligner(seq1, seq2, matchScore, mismatchScore, gapopen, gapext);
         aligner.setIdentityCutoff(0.7);
         numnodes.push_back(aligner.run());
         if (aligner.getIdentity()>maxidentity) {
            maxidentity = aligner.getIdentity();
            vector<int> tmp={gapopen, gapext, numnodes.back()};
            result.clear();
            result.insert(tmp);
         }
         else if (aligner.getIdentity() == maxidentity) {
            vector<int> tmp={gapopen, gapext, numnodes.back()};
            result.insert(tmp);
         }
      }
   }
   time(&timer2);
   double seconds = difftime(timer2, timer1);
   stddev stat;
   for (int i=0; i<numnodes.size(); ++i) {
      stat(numnodes[i]);
   }
   vector<string> header={"match", "mismatch", "gapopen", "gapext", 
      "nnExplored", "identity", "time", "numTest", 
      "nnE_min", "nnE_max", "nnE_mean", "nnE_std"};
   copy(header.begin(), header.end(), ostream_iterator<string>(cout, "\t"));
   cout << endl;
   //cout << result.size() << " best result sets\n";
   for (auto ptr=result.begin(); ptr != result.end(); ++ptr) {
      cout << matchScore << '\t' << mismatchScore << '\t';
      for (int j=0; j<(*ptr).size(); ++j) {
         cout << (*ptr)[j] << '\t';
      }
      cout << maxidentity << '\t' << seconds << '\t' << trial1*trial2
         << '\t' << *min_element(numnodes.begin(), numnodes.end()) 
         << '\t' << *max_element(numnodes.begin(), numnodes.end())
         << '\t' << stat.getMean() << '\t' << stat.getStd() << endl;
   }
   //cout << "number of nodes explored in each round:\n";
   //copy(numnodes.begin(), numnodes.end(), ostream_iterator<int>(cout, ", "));
   //cout << string(50, '=') << endl;
}
