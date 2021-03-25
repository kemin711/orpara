#include "bioseq.h"
#include "scorematrix.h"
#include "dynalnt.h"
#include "fastq.h"
#include <strformat.h>
#include <iostream>
#include <cstring>
#include <string>

/** a program to tes the alignemnt classes.
 */
using namespace std;
using namespace orpara;

class ProgParam {
   private:
      float identitycut;
      int alnlencut;
      /**
       * Reverse complement 1 for the reference. 
       * 2 for database, 3 for plus strand of database,
       * and 0 for try both strand of query database.
       */
      int revcomp;

   public:
      ProgParam() 
         : identitycut(0.75), alnlencut(30), revcomp(0) { }
      void setAlnlenCut(const int al) {
         alnlencut = al;
      }
      void setIdentityCut(const float ic) {
         identitycut = ic;
      }
      float getIdentityCut() const {
         return identitycut;
      }
      int getAlnlenCut() const {
         return alnlencut;
      }
      void reverseComplementReference() {
         revcomp = 1;
      }
      /**
       * Costly 
       */
      void reverseComplementDatabase() {
         revcomp = 2;
      }
      void setRevcomp(int rv) {
         revcomp=rv;
      }
      int useReverseComplement() const {
         return revcomp;
      }
      bool useBothStrand() const {
         return revcomp == 0;
      }
      bool usePlusStrand() const {
         return revcom == 3;
      }
      void setPlusStrand() {
         revcom = 3;
      }
      int getRevcomp() const {
         return revcomp;
      }
};

void usage() {
   cerr << "alnlocalmany -r REFERENCE -o OUTPUT DBFILE\n"
      << "Options:\n"
      << "   -r reference sequence file only one sequence allowed\n"
      << "   -o output file name if not given this program will make\n"
      << "       one for you\n"
      << "   -d database file contain many sequences in fasta or fastq format\n"
      << "       The program will judge the type of sequence based on the \n"
      << "       file suffix. *.fastq files will be treated as fastq files\n"
      << "       and others will be treated as fasta files\n"
      << "   --identity-cut FLOAT a fraction number (0, 1] to filter\n"
      << "     alignments whose identity < this will be discarded\n"
      << "   --alnlen-cut INTEGER a integer value (default 30) to filter\n"
      << "     alignment. If alignment length < this value, then will be discarded\n"
      << "   --minus-reference flag to use negative strand of reference\n"
      << "   --minus-database flag to use negative strand of database\n"
      << "     this is a costly operation, every sequence in db will be\n"
      << "     reverse complemented before alignment\n"
      << "   --plus-database use plus strand of data file. Default use both\n"
      << "   --help will print this message\n"
      << "Positional argument:\n"
      << "   the database file can be given as the only positional argument\n"
      << "   if it is in fasta foramt. If input in fastq format, you must\n"
      << "   use the -q option\n";
}

string makeOutputFile(const string& ref, const string& db);
bool isDNA(const string& infile);
bool isProtein(const string& infile);
void alignDNAMany(const string& ref, const string& dbf, const string& outf, const ProgParam& param);
void alignDNAManyFastq(Dynaln<SimpleScoreMethod>& matcher, 
      ifstream& inf, ofstream& ouf, const ProgParam &param, int& sqcnt);
void alignProteinMany(const string& ref, const string& dbf, const string& outf);

/**
 * Aling many sequences with a single reference sequence 
 * and generate a summary file.
 */
int main(int argc, char *argv[]) {
   string reffile, dbfile, outfile;
   //reffile="chr3_178935841_178936341.fas";
   //dbfile="CF_HD786_Old_3_prey_R1F.fastq";
   ProgParam param;
   int i=1;
   while (i < argc) {
      if (!strcmp(argv[i], "--help")) {
         usage();
         return 0;
      }
      else if (!strcmp(argv[i], "-r")) {
         reffile=argv[++i];
      }
      else if (!strcmp(argv[i], "-d")) {
         dbfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         outfile = argv[++i];
      }
      //else if (!strcmp(argv[i], "-q")) { dbfile=argv[++i]; }
      else if (!strcmp(argv[i], "--identity-cut")) { 
         param.setIdentityCut(atof(argv[++i])); 
      }
      else if (!strcmp(argv[i], "--alnlen-cut")) { 
         param.setAlnlenCut(atoi(argv[++i])); 
      }
      else if (!strcmp(argv[i], "--minus-reference")) { 
         param.reverseComplementReference(); 
      }
      else if (!strcmp(argv[i], "--minus-database")) { 
         param.reverseComplementDatabase(); 
      }
      else if (!strcmp(argv[i], "--plus-database")) {
         param.setPlusStrand();
      }
      else {
         dbfile = argv[i];
      }
      ++i;
   }
   if (reffile.empty()) {
      usage();
      cerr << "reference file not provided\n";
      return 1;
   }
   if (dbfile.empty()) {
      usage();
      cerr << "must provide a data file of one or more sequences in fasta or fastq format\n";
      return 1;
   }
   if (outfile.empty()) {
      outfile=makeOutputFile(reffile, dbfile);
   }
   if (isDNA(reffile)) {
      alignDNAMany(reffile, dbfile, outfile, param);
   }
   else if (isProtein(reffile)) {
      alignProteinMany(reffile, dbfile, outfile);
   }
   else {
      cerr << "reference of unknown type, failed\n";
      return 1;
   }

   return 0;
}

bool isDNA(const string& infile) {
   bioseq bsq;
   if (!bsq.read(infile)) {
      throw runtime_error("failed to get reference seq");
   }
   if (bsq.guessType() == DNASEQ) {
      return true;
   }
   return false;
}
bool isProtein(const string& infile) {
   bioseq bsq;
   if (!bsq.read(infile)) {
      throw runtime_error("failed to get reference seq");
   }
   if (bsq.guessType() == PROTEINSEQ) {
      return true;
   }
   return false;
}

// TODO: needs to be tested
void alignDNAMany(const string& ref, const string& dbf, const string& outf, const ProgParam& param)
{
   // align ref to fasta file
   SimpleScoreMethod sm(10, -9, -29, -3);
   Dynaln<SimpleScoreMethod> aligner(sm), revaligner(sm);
   DNA dnaref;
   dnaref.read(ref);
   if (param.useReverseComplement() == 1) {
      dnaref.revcomp();
   }
   aligner.setSeq1(dnaref);
   revaligner.setSeq1(dnaref);
   ifstream inf(dbf);
   if (inf.fail()) {
      throw runtime_error("failed to open dbfile " + dbf);
   }
   ofstream ouf(outf);
   if (ouf.fail()) {
      throw runtime_error("failed to open outpufle " + outf);
   }
   Dynaln<SimpleScoreMethod>::printSummaryHeader(ouf, "\t", false);
   ouf << endl;

   int numseq=0;
   // diverge according to different methods
   if (endwith(dbf, ".fastq")) {
      alignDNAManyFastq(aligner, inf, ouf, param, numseq);
   }
   else {
      DNA dna; 
      while (dna.read(inf)) {
         if (param.useBothStrand()) { // == 0
            aligner.setSeq2(dna);
            DNA rcdna = dna.revcompCopy();
            revaligner.setSeq2(rcdna);
            int score1 = aligner.loca();
            int score2 = revaligner.local();
            if (score1 >= score2) {
               if (aligner.getIdentity() > param.getIdentityCut() 
                     && aligner.getAlnlen() > param.getAlnlenCut()) 
               {
                  aligner.buildResult(0,0);
                  aligner.printSummary(ouf, "\t", false) << endl;
               }
            }
            else {
               if (revaligner.getIdentity() > param.getIdentityCut()
                     && revaligner.getAlnlen() > param.getAlnlenCut()) {
                  revaligner.buildResult(0,0);
                  revaligner.printSummary(ouf, "\t", false) << endl;
               }
            }
         }
         else { // puls or minus
            if (param.useReverseComplement() == 2) {
               dna.revcomp();
            }
            aligner.setSeq2(dna);
            aligner.runlocal();
            if (aligner.getIdentity() > param.getIdentityCut() 
                  && aligner.getAlnlen() > param.getAlnlenCut()) 
            {
               aligner.printSummary(ouf, "\t", false) << endl;
            }
         }
         ++numseq;
      }
   }
   cerr << numseq << " database sequences aligned to reference, result in "
      << outf << endl;
}

void alignDNAManyFastq(Dynaln<SimpleScoreMethod>& matcher, ifstream& inf, ofstream& ouf, const ProgParam &param, int& sqcnt) {
   Fastq fsq;
   while (fsq.read(inf)) {
      DNA raw(fsq.getName(), fsq.getSequence());
      if (param.useReverseComplement() == 2) {
         raw.revcomp();
      }
      matcher.setSeq2(raw);
      matcher.runlocal();
      if (matcher.getIdentity() > param.getIdentityCut() 
            && matcher.getAlnlen() > param.getAlnlenCut()) 
      {
         matcher.printSummary(ouf, "\t", false) << endl;
      }
      ++sqcnt;
   }
}

void alignProteinMany(const string& ref, const string& dbf, const string& outf) {
   Dynaln<ProteinScoreMethod> alnp;
   Protein pref;
   if (!pref.read(ref)) {
      throw runtime_error("Failed to read protein reference file: " + ref);
   }
   alnp.setSeq1(pref);
   //alnp.setGapParameter(gapOpen, gapExtend);
   ifstream inf(dbf);
   if (inf.fail()) {
      throw runtime_error("Failed to open protein database file " + dbf);
   }
   ofstream ouf(outf);
   if (ouf.fail()) {
      throw runtime_error("Failed to open output file " + outf);
   }
   Protein prt;
   while (prt.read(inf)) {
      alnp.setSeq2(prt);
      alnp.runlocal(); // seqXbegin for labeling the alignment
      //alnp.printAlign(ous);
      alnp.printSummary(ouf, "\t", false) << endl;
   }
   cerr << "protein alignment written to " << outf << endl;
}

string makeOutputFile(const string& ref, const string& db) {
   string out = getFileStem(ref);
   //string out = "blala";
   out += "_";
   out += getFileStem(db);
   out += ".aln";
   return out;
}
