#include "matrix.h"
#include "bioseq.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
/* No longer use this one.
 * Migrated to scorematrix.h
 */

using namespace std;


int main(int argc, char *argv[]) {
   Matrix pm("blosum62.50");
   int neighborThreshold=10;
   float neighborFraction=0.75;
   int i=1;

   while (i<argc) {
      if (!strcmp(argv[i], "-f")) neighborFraction=atof(argv[++i]);
      else if (!strcmp(argv[i], "-t")) neighborThreshold=atoi(argv[++i]);
      else {
         cerr << argv[i] << " is not recognized\n";
      }
      ++i;
   }

   /*
   vector<string> wds;
   pm.getWords(wds, 3);
   for (int i=0; i<wds.size(); i++) {
      cout << wds[i] << " ";
      if ((i+1)%10 == 0) cout << endl;
   }
   cout << endl;
   cout << "vector<string> version tested fine\n";
   */

   char **aw=pm.allwords(2);
   ofstream outf("neighbors.txt");
   ofstream BAD("noneighbor.txt");
   //pm.showWords();
   map<int, int> stat;
   int numn;
   for (i=0; i<pm.getNumberOfWords(); i++) {
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
   Matrix default_pm;
   default_pm.show();

   Matrix nuc;
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
   nuc.setMatrix("NUC.4.4");
   cout << "min score: " << nuc.getMinScore() << endl;
   cout << "max score: " << nuc.getMaxScore() << endl;
   nuc.show();
   
   return 0;
}
