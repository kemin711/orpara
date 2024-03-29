#include <iostream>
#include "stddev.h"
#include <vector>
#include <cstring>

using namespace std;
using namespace orpara;

/**
 * For testing partition function from stddev.h
 */
int main(int argc, char* argv[]) {
   if (argc == 1) {
      cerr << "testing no need for any arguments\n";
   }
   int i=1;
   while (i < argc) {
      if (!strcmp(argv[i], "--help")) {
         cerr << "just a test\n";
      }
      ++i;
   }
   vector<unsigned long> input{476232019, 487787336, 490498321, 500103870, 529907050, 540695814, 541056067, 546623159,
         552267810, 552818236, 557372854, 560951249, 565116326, 568780772, 574660216, 578322332,
         582110850, 595405819, 604857617, 618781588, 621258109, 623640992, 632250899, 632329067,
         636619716, 637425233, 639571346, 646165702, 650282086, 652572857, 652618365, 657486829,
         657637713, 664078742, 672483244, 673244078, 680649430, 692453644, 694445694, 694541912,
         708384169, 710398793, 720181320, 735814853, 15029402817, 15368636177, 27564355721, 28210994308};
   Bisectsorted<unsigned long> cutter(std::move(input));
   cutter.separate();
   cout << endl << "H/L: " << cutter.getHighLowRatio() << " log10(H/L): " 
      << cutter.getLog10Ratio() << " pivot: " << cutter.getPivot() << endl
      << "total ratio: " << cutter.getTotalHighLowRatio() 
      << "total log ratio: " << cutter.getTotalLog10Ratio() << endl;
}

