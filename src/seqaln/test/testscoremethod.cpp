#include "scorematrix.h"
#include "bioseq.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

void testScoreMethod() {
   ScoreMethod sm;
   cout << "match: " << sm.getMatch() << " mismatch: " << sm.getMismatch() << endl;
   sm.show(cout);
}

void testProteinMethod() {
   cout << "Testing proteinScoreMethod ...\n";

   ProteinScoreMethod pm("blosum62.50");
   cout << "created ProteinScoreMethod instance\n";
   pm.show(cerr);
   
   vector<string> wds;
   pm.getWords(wds, 3);
   for (int i=0; i<wds.size(); i++) {
      cout << wds[i] << " ";
      if ((i+1)%10 == 0) cout << endl;
   }
   cout << endl;
   cout << "vector<string> version tested fine\n";

   char **aw=pm.allwords(2);
   ofstream outf("neighbors.txt");
   ofstream BAD("noneighbor.txt");
   int neighborThreshold=10;
   float neighborFraction=0.75;
   //pm.showWords();
   map<int, int> stat;
   int numn;
   for (int i=0; i<pm.getNumberOfWords(); i++) {
      numn = pm.similarWord_debug(outf, aw[i], 2, neighborThreshold, neighborFraction);
      ++(stat[numn]);
      if (numn == 0) {
         BAD << aw[i] << endl;
      }
   }
   map<int,int>::const_iterator it;
   cout << "#number of neighbor | count at threshold=" << neighborThreshold
      << " fraction=" << neighborFraction << "\n";
   for (it=stat.begin(); it != stat.end(); it++) {
      cout << it->first << "\t" << it->second << endl;
   }

   cerr << "Testing default matrix\n";
   ProteinScoreMethod default_pm;
   default_pm.show(cerr);
   cerr << default_pm.lookup('S', 'T') << " score for S x T\n";
   cerr << default_pm.lookup('Y', 'Y') << " score for Y x Y\n";
   cerr << "score for gap open: " << default_pm.getGapOpen() << endl;
   cerr << "score for gap extend: " << default_pm.getGapExtend() << endl;
}

void testNucMethod() {
   NucleicScoreMethod nuc;
   char* homedir = getenv("HOME");
   if (homedir != NULL) {
      cout << "My home: " << homedir << endl;
   }
   else {
      cout << "Failed to get home dir\n";
   }
   string matrixPath(homedir);
   matrixPath += "/src/proj/seqaln/matrix";
   //nuc.setPath("/home/zhouke/scr/proj/seqaln/matrix
   nuc.setPath(matrixPath);
   //nuc.setMatrix("NUC.4.4");
   nuc.use("NUC.4.4");
   cout << "min score: " << nuc.getMinScore() << endl;
   cout << "max score: " << nuc.getMaxScore() << endl;
   nuc.show();
}

template<class T>
void testAny(T &sm) {
   cout << "AxA from NucleicScoreMethod: " << sm.lookup('A', 'A')
      << "\nNxN from ProteinScoreMethod: " << sm.lookup('N', 'N')
      << "\nWxW from ProteinScoreMethod: " << sm.lookup('W', 'W')
      << "\nSxS from ProteinScoreMethod: " << sm.lookup('S', 'S') << endl;
}

void testMixed() {
   cout << "Testing mixed operations ...\n";
   ScoreMethod sm;
   ProteinScoreMethod psm;
   psm.show();
   NucleicScoreMethod nsm;
   cout << "AxA from ScoreMethod: " << sm.lookup('A', 'A') 
      << "\nAxA from NucleicScoreMethod: " << nsm.lookup('A', 'A')
      << "\nNxN from ProteinScoreMethod: " << psm.lookup('N', 'N')
      << "\nSxS from ProteinScoreMethod: " << psm.lookup('S', 'S') << endl;
   cout << "N's code: " << aachar2num('N') << endl;
   cout << "S's code: " << aachar2num('S') << endl;

   ScoreMethod *p = new ProteinScoreMethod;
   cout << "WxW from ScoreMethod pointer to protein: " << p->lookup('W', 'W') << endl;
   delete p;

   // testing the template method
   testAny(sm);
   testAny(psm);
}



int main(int argc, char *argv[]) {
   int i=1;
   int neighborThreshold=10;
   float neighborFraction=0.75;

   while (i<argc) {
      if (!strcmp(argv[i], "-f")) neighborFraction=atof(argv[++i]);
      else if (!strcmp(argv[i], "-t")) neighborThreshold=atoi(argv[++i]);
      else {
         cerr << argv[i] << " is not recognized\n";
      }
      ++i;
   }

   testScoreMethod();
   testProteinMethod();
   testNucMethod();
   testMixed();
   
   return 0;
}


