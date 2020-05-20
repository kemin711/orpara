#include <iostream>
#include "dynalnt.h"


using namespace std;
using namespace orpara;

int main(int argc, char* argv[]) {
   string sA1 = "ACCCGTGTGGCCAGTCGTACGTGTACACTGACGTACGATCAGTC";
   string sB1 = "ACGTGTCCAGTCCGTGCTGACGAATTGTGTGTTTTGGTTTC";
   string sA2 = "ACGGGTGTGCCCAGACGT";
   string sB2 = "TTTTTTGGGGTTTGGC";

   SimpleScoreMethod sm(5, -4, -8, -8);
   int maxLen = max(max(sA1.size(), sA2.size()), max(sB1.size(), sB2.size()));

   cout << "seqA : " << sA1 << " -- " << sA2 << endl;
   cout << "seqB : " << sB1 << " -- " << sB2 << endl;

   Dynaln<SimpleScoreMethod> aln(sm);
   double score12, score11, score22;
   vector<double> normalScore(maxLen+1);

   for (int l = 1; l <= maxLen; l++) {
   //for (int l = 3; l <= maxLen; l++) {
      int x=max<int>(0, int(sA1.size())-l);
      int y=max<int>(0, int(sB1.size()) - l);
      cout << string(2 * maxLen + 20, '=') << endl;
      cout << "\nl=" << l << " x y: " << x << " " << y << endl;

      DNA parSeqA1("seqA1", sA1.substr(x, l));
      DNA parSeqB1("seqB1", sB1.substr(y, l));
      DNA parSeqA2("seqA2", sA2.substr(0, l));
      DNA parSeqB2("seqB2", sB2.substr(0, l));

      aln.setSeq(parSeqA1, parSeqB1);
      aln.runlocal(); // if only interested in score, then use local() to save computation
      score11 = aln.getScore();
      cout << parSeqA1 << parSeqB1 << endl;
      aln.printAlign(cout);
      cout << string(30, '&') << endl;

      aln.setSeq2(parSeqB2);
      aln.runlocal();
      score12 = aln.getScore();
      cout << parSeqA1 << parSeqB2 << endl;
      aln.printAlign(cout);
      cout << string(30, '&') << endl;

      aln.setSeq1(parSeqA2);
      aln.runlocal();
      score22 = aln.getScore();
      cout << parSeqA2 << parSeqB2 << endl;
      aln.printAlign(cout);
      cout << string(30, '&') << endl;
      cout << "score 1x1=" << score11 << " 1x2=" << score12 << " 2x2=" << score22 << endl;
   }
   return 0;
}

