#include "scorematrix.h"
#include "bioseq.h"
#include "dynalnt.h"
#include <fstream>
#include "alnexaminer.h"

void runTest2(const string &ref, const string &read);
void runTest();

void extractSegment(const Dynaln<SimpleScoreMethod> &aln, 
      const Alnexaminer &zoner, const DNA &seq) {
   vector<Alnseg> zone=zoner.getSegment();
   vector<pair<int,int> > matchIndex=aln.getAlnindexVector(); 
   for (size_t i=0; i<zone.size(); ++i) {
      if (zone[i].identity > 0.99) {
         // begin index of second sequence
         int b=matchIndex[zone[i].b].second;
         int e=matchIndex[zone[i].e-1].second;
         cout << b << "-" << e << endl;
         cout << seq.subseq(b+1, e+1) << endl;
      }
      else {
         cout << zone[i] << endl
            << " not good enough\n";
      }
   }
}

/** 
 * a program to test the alignemnt classes.
 */
int main(int argc, char *argv[]) {
   string refseq, read;
   int i=1;
   while (i < argc) {
      if (!strcmp(argv[i], "-f")) {
         refseq=argv[++i];
      }
      else if (!strcmp(argv[i], "-r")) {
         read=argv[++i];
      }
      else {
         if (i+1 < argc) {
            refseq=argv[i];
            read=argv[i+1];
            ++i;
         }
         else {
            cerr << "bad arguments\n";
            exit(1);
         }
      }
      ++i;
   }
   runTest2(refseq, read);
   return 0;
}

// test only two files
void runTest2(const string &ref, const string &read) {
   cerr << "testing from " << ref << " and " << read << endl;
   DNA dnaRef, dnaRead;
   dnaRef.read(ref);
   dnaRead.read(read);
   SimpleScoreMethod spcm(13, -11, 29, 3);
   Dynaln<SimpleScoreMethod> aln(spcm);
   aln.setSeq(dnaRef, dnaRead);
   aln.runlocal();
   aln.printAlign(cout);
   Alnexaminer zoner;
   zoner(aln.getMiddleAln());
   zoner.debugShow(cout);
   if (zoner.isChimera()) {
      zoner.printResult(cout) << endl;
      extractSegment(aln, zoner, dnaRead);
   }
   else {
      cout << "is not chimera\n";
   }
}

void runTest() {
   cerr << "Testing chimera detection ...\n";
   string referenceFile("cluster0.fas");
   string chimeraFile("chimera.fas");
   SimpleScoreMethod spcm(13, -11, 29, 3);
   Dynaln<SimpleScoreMethod> aln(spcm);
   DNA ref, chimera;
   ref.read(referenceFile);
   ifstream inf(chimeraFile);
   if (inf.fail()) {
      throw runtime_error("cannot locate file: " + chimeraFile);
   }
   aln.setSeq1(ref);
   string line;
   Alnexaminer zoner;
   DNAQualCount* sgl;
   while (chimera.read(inf)) {
      aln.setSeq2(chimera);
      aln.runlocal();
      aln.printAlign(cout);
      zoner(aln.getMiddleAln());
      if (zoner.isChimera()) {
         zoner.printResult(cout) << endl;
         extractSegment(aln, zoner, chimera);
      }
      else {
         cout << "is not chimera\n";
      }
      cout << "alignment length: " << aln.getAlnlen() << endl;
   }
}



