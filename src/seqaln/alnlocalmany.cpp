#include "bioseq.h"
#include "scorematrix.h"
#include "dynalnt.h"
#include "fastq.h"


/** a program to tes the alignemnt classes.
 */
using namespace std;
using namespace orpara;

void usage() {
   cout << "alnlocal seq1.fas seq2.fastq -o outfile\n"
      << "Options:\n"
      << "        -r [1,2]  reverse complement first or second sequence\n"
      << "        -f reference sequence file\n"
      << "        -q fastq sequence file containing many input sequences\n"
      << "        -b1 1-based index of sequence1 begin\n"
      << "        -e1 1-based index of sequence1 end\n"
      << "        -b2 1-based index of sequence2 begin\n"
      << "        -e2 1-based index of sequence2 end\n"
      << "        --gap-open or -g  gap open cost, default min(matrix score)*2 \n"
      << "        --gap-extend or -e  gap extend cost, default min(matrix score)\n"
      << "        -o output file, in summary tabular format\n"
      << "         if output file name is not specified, the program will generate one\n";
}

/**
 * Remove the path part
 * Only works for unix path separator '/'.
 */
string getFilePart(const string &path);

/**
 * Remove the path part and suffix
 * @return the file stem only.
 */
string getOnlyFileStem(const string &fname);
/**
 * @param b begin of the sequence
 * @param e end of the sequence
 * @seq will be altered to subsequence of [b, e] in 1-based index.
 */
void shrinkSequence(bioseq &seq, int b, int e);
//void setGapParameters(Matrix &scm, int go, int ge);
template<class T>
void showAlignInfo(const Dynaln<T> &aln, ostream &ous);
void alignProtein(const int gapOpen, const int gapExtend, const bioseq &p1, const bioseq &p2, 
      const int seq1begin, const int seq2begin, ostream &ous);
void alignDNA(const int gapo, const int gape, const DNA &dna1, const DNA &dna2,
      const int seq1begin, const int seq2begin, ostream &ous);
void alignSimple(const bioseq &s1, const bioseq &s2,
      const int seq1begin, const int seq2begin, ostream &ous);

int main(int argc, char *argv[]) {
   int i = 1;
   string file1, file2, outfile;
   int reverseComplement = 0;
   int seq1begin=1, seq1end=-1, seq2begin=1, seq2end=-1;
   int gapOpen = 1;
   int gapExtend = 1;
   bool simpleMethod = true;
   while (i < argc) {
      if (string(argv[i]) == "-f") { file1=argv[++i]; }
      else if (string(argv[i]) == "-q") { file2=argv[++i]; }
      else if (string(argv[i]) == "-o") { outfile=argv[++i]; }
      else if (string(argv[i]) == "-r") { reverseComplement=atoi(argv[++i]); }
      else if (string(argv[i]) == "-b1") { seq1begin=atoi(argv[++i]); }
      else if (string(argv[i]) == "-e1") { seq1end=atoi(argv[++i]); }
      else if (string(argv[i]) == "-b2") { seq2begin=atoi(argv[++i]); }
      else if (string(argv[i]) == "-e2") { seq2end=atoi(argv[++i]); }
      else if (string(argv[i]) == "--gap-open" 
            || string(argv[i]) == "-g") { gapOpen=atoi(argv[++i]); }
      else if (string(argv[i]) == "--gap-extend"
            || string(argv[i]) == "-e") { gapExtend=atoi(argv[++i]); }
      else {
         file1 = argv[i];
         if (i+1 < argc && argv[i+1][0] != '-') {
            file2 = argv[++i];
         }
      }
      ++i;
   }

   if (file1.empty() || file2.empty()) {
      usage();
      return 1;
   }
   if (outfile.empty() && !file1.empty() && !file2.empty()) {
      outfile = getOnlyFileStem(file1);
      outfile += "_";
      outfile += getOnlyFileStem(file2);
      if (seq1begin > 1 || seq2begin > 1 || seq1end != -1 || seq2end != -1) {
         outfile += "sub";
      }
      outfile += ".aln";
   }
   ofstream ouf(outfile.c_str());
   if (!ouf) {
      cerr << "Failed to open " << outfile << " for alignment output!\n";
      return 1;
   }

   bioseq seq1, seq2;
   seq1.read(file1); // reference
   //seq2.read(file2);
   ifstream inf(file2.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << file2 << endl;
      return 1;
   }

   Fastq fastq;

   if (seq1.guessType() == DNASEQ) {
      SimpleScoreMethod sm(10, -9, -29, -3);
      Dynaln<SimpleScoreMethod> aligner(sm);
      DNA dnaref = DNA(seq1);
      aligner.setSeq1(dnaref);

      Dynaln<SimpleScoreMethod>::printSummaryHeader(ouf, "\t", false);
      ouf << endl;
      while (fastq.read(inf)) {
         DNA raw(fastq.getName().substr(1), fastq.getSequence());
         aligner.setSeq2(raw);
         aligner.runlocal();
         aligner.printSummary(ouf, "\t", false) << endl;
      }
   }
   else { // protein align
      if (simpleMethod) 
         alignSimple(seq1, seq2, seq1begin, seq2begin, ouf);
      else 
         alignProtein(gapOpen, gapExtend, seq1, seq2, seq1begin, seq2begin, ouf);
   }
   cout << "aligment written to file: " << outfile << endl;

   return 0;
}

void alignSimple(const bioseq &s1, const bioseq &s2,
      const int seq1begin, const int seq2begin, ostream &ous) {
   SimpleScoreMethod sm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethod> aln(s1, s2);
   aln.setMatrix(sm);
   aln.runlocal(seq1begin-1, seq2begin-1);
   aln.printAlign(ous);
   showAlignInfo<SimpleScoreMethod>(aln, cout);
}

void alignProtein(const int gapOpen, const int gapExtend, 
      const bioseq &p1, const bioseq &p2, 
      const int seq1begin, const int seq2begin, 
      ostream &ous) 
{
   // protein align using default matrix
   Dynaln<ProteinScoreMethod> alnp(p1, p2);
   if (gapOpen < 1) 
      alnp.setGapParameter(gapOpen, gapExtend);
   alnp.runlocal(seq1begin-1, seq2begin-1); // seqXbegin for labeling the alignment
   alnp.printAlign(ous);
   showAlignInfo<ProteinScoreMethod>(alnp, cout);
}

void alignDNA(const int gapo, const int gape, const DNA &dna1, const DNA &dna2,
      const int seq1begin, const int seq2begin, ostream &ous) {
   Dynaln<NucleicScoreMethod> aln(dna1, dna2);
   if (gapo < 1 || gape < 1) {
      aln.setGapParameter(gapo, gape);
   }
   aln.runlocal(seq1begin-1, seq2begin-1);
   aln.printAlign(ous);
   showAlignInfo<NucleicScoreMethod>(aln, cout);
}

template<class T>
void showAlignInfo(const Dynaln<T> &aln, ostream &ous) {
   ous << aln.getSeq1Name() << "|" << aln.getSeq2Name()
      << " seqlength=" << aln.getSeq1Length() << "," << aln.getSeq2Length()
      << " alnlength=" << aln.getAlnlen() 
      << " numgap=" << aln.getNumgaps1() << "," << aln.getNumgaps2()
      << " gaplength=" << aln.getGaplen1() << "," << aln.getGaplen2()
      << " identity=" << aln.getIdentity() 
      << " range1=" << aln.topBeginIndex()+1 << "-" << aln.topEndIndex()+1
      << " range2=" << aln.bottomBeginIndex()+1 << "-" << aln.bottomEndIndex()+1
      << endl;
}

string getFilePart(const string &path) {
   string::size_type i = path.rfind('/');
   if (i != string::npos) {
      ++i;
      return path.substr(i);
   }
   else return path;
}

string getOnlyFileStem(const string &fname) {
   string filepart = getFilePart(fname);
   return getFileStem(filepart);
}


void shrinkSequence(bioseq &seq, int b, int e) {
   if (b > 1) {
      seq = seq.subseq(b, e);
   }
   else {
      if (e != -1) {
         seq = seq.subseq(1, e);
      }
   }
}

/*
void setGapParameters(Matrix &scm, int go, int ge) {
   if (go != 0) {
      if (go > 0) go = -go;
      scm.setGapInsert(go);
   }
   if (ge != 0) {
      if (ge > 0) ge = -ge;
      scm.setGapExtend(ge);
   }
}
*/
