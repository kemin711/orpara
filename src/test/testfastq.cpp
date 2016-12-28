#include "fastq.h"
#include <string>
#include <iostream>
#include <fstream>


using namespace std;

/**
 * testfastq test.fastq
 */
int main(int argc, char* argv[]) {
   string infile="test.fastq";
   if (argc > 1) infile=argv[1];
   ifstream inf(infile);
   if (inf.fail()) {
      cerr << "cannot open " << infile << endl;
      return 1;
   }

   Fastq fq;
   int cnt=0;
   while (fq.read(inf)) {
      ++cnt;
      cout << fq.getName() << endl;
   }
   cout << cnt << " sequences in the file\n";

   return 0;
}

