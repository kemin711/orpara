#include <iostream>
#include "scorematrix.h"
#include "bioseq.h"
#include "dynalnt.h"

/** a program to tes the alignemnt classes.
 */
using namespace std;
using namespace orpara;

/*
template<class T>
void showAlignInfo(const Dynaln<T> &aln, ostream &ous) {
   ous << aln.getSeq1Name() << "|" << aln.getSeq2Name()
      << endl << " seqlength=" << aln.getSeq1Length() << "," << aln.getSeq2Length()
      << endl << " alnlength=" << aln.getAlnlen() 
      << endl << " numgap=" << aln.getNumgaps1() << "," << aln.getNumgaps2()
      << endl << " gaplength=" << aln.getGaplen1() << "," << aln.getGaplen2()
      << endl << " identity=" << aln.getIdentity() 
      << endl << " range1=" << aln.topBeginIndex()+1 << "-" << aln.topEndIndex()+1
      << endl << " range2=" << aln.bottomBeginIndex()+1 << "-" << aln.bottomEndIndex()+1
      << endl;
}

void plotScores(double* normalScore, int maxLen) {
   string bar = string(maxLen, '.');
   string spc = string(maxLen, ' ');

   vector<string> graph;

   for (int i = 0; i <= 50; i++)
       if (i % 10 == 0)
           graph.push_back(bar);
       else
           graph.push_back(spc); 
   
   for (int l = 1; l <= maxLen; l++)
       graph[10 * normalScore[l]][l] = '*';
   
   for (int i = 0; i <= 50; i++)
       cout << graph[50 - i] << endl;
}


void align(const string &sA1, const string &sA2, const string &sB1, const string &sB2, 
      int match, int mismatch, int gapOpen, int gapExtend) 
{
   SimpleScoreMethod sm(match, mismatch, gapOpen, gapExtend);
   int maxLen = max(max(sA1.size(), sA2.size()), max(sB1.size(), sB2.size()));

   cout << "seqA : " << sA1 << " -- " << sA2 << endl;
   cout << "seqB : " << sB1 << " -- " << sB2 << endl;

   Dynaln<SimpleScoreMethod> aln(sm);
   double normalScore[maxLen + 1], score12, score11, score22;

   //for (int l = 1; l <= maxLen; l++) {
   for (int l = 3; l <= maxLen; l++) {
      int x, y;
      x=max<int>(0, int(sA1.size())-l);
      y=max<int>(0, int(sB1.size()) - l);
      if (x > sA1.size()-1 || y > sB1.size()) {
         throw runtime_error("x, y passed end of input string, cannot take substring beond end");
      }
      cout << string(2 * maxLen + 20, '=') << endl;
      cout << "\nl=" << l << " x y: " << x << " " << y << endl;

      DNA parSeqA1("seqA1", sA1.substr(x, l));
      DNA parSeqB1("seqB1", sB1.substr(y, l));
      DNA parSeqA2("seqA2", sA2.substr(0, l));
      DNA parSeqB2("seqB2", sB2.substr(0, l));

      aln.setSeq(parSeqA1, parSeqB1);
      aln.runlocal();
      score11 = aln.getScore();
      cout << parSeqA1 << parSeqB1 << endl;
      aln.printAlign(cout);
      cout << string(30, '&') << endl;

      //cout << "parSeqA1: " << parSeqA1.toString() << endl;
      //cout << "parSeqB2: " << parSeqB2.toString() << endl;
      aln.setSeq2(parSeqB2);
      aln.runlocal();
      score12 = aln.getScore();
      cout << parSeqA1 << parSeqB2 << endl;
      aln.printAlign(cout);
      cout << string(30, '&') << endl;

      //cout << "parSeqA2: " << parSeqA2.toString() << endl;
      //cout << "parSeqB2: " << parSeqB2.toString() << endl;
      aln.setSeq1(parSeqA2);
      aln.runlocal();
      score22 = aln.getScore();
      cout << parSeqA2 << parSeqB2 << endl;
      aln.printAlign(cout);
      cout << string(30, '&') << endl;

      cout << "score 1x1=" << score11 << " 1x2=" << score12 << " 2x2=" << score22 << endl;
      //normalScore[l] = max(max(score12, score22), score11) /  static_cast<double>(l);
      //cout << normalScore[l];
   }

   //aln.printAlign(cout);
   //showAlignInfo<SimpleScoreMethod>(aln, cout);
   
   //plotScores(normalScore, maxLen);
}

*/

int main(int argc, char *argv[])
{
   string sA1 = "ACCCGTGTGGCCAGTCGTACGTGTACACTGACGTACGATCAGTC";
   string sB1 = "ACGTGTCCAGTCCGTGCTGACGAATTGTGTGTTTTGGTTTC";
   string sA2 = "ACGGGTGTGCCCAGACGT";
   string sB2 = "TTTTTTGGGGTTTGGC";

   int match = 5, mismatch = -4, gapOpen = -8 , gapExtend = -8; 

   //align(sA1, sA2, sB1, sB2, match, mismatch, gapOpen, gapExtend);
    
   return 0;
}

