#include "scorematrix.h"
#include "bioseq.h"
#include "dynalnt.h"
#include <ctime>

/** a program to tes the alignemnt classes.
 */
using namespace std;
using namespace orpara;

struct Progparam {
   public:
      Progparam() : m(5), s(-3), g(-23), e(-1) {}
      Progparam(int match, int mismatch, int gapo, int gape) 
         : m(match), s(mismatch), g(gapo), e(gape) {}

      int m,s,g,e;
};


void usage() {
   cout << "alnglobal seq1.file seq2.file -o outfile\n"
      << "Options:\n"
      << "        -r [1,2]  reverse complement first or second sequence\n"
      << "        -1 first sequence file\n"
      << "        -2 second sequence file\n"
      << "        -b1 1-based index of sequence1 begin\n"
      << "        -e1 1-based index of sequence1 end\n"
      << "        -b2 1-based index of sequence2 begin\n"
      << "        -e2 1-based index of sequence2 end\n"
      << "        --gap-open or -g  gap open cost, default min(matrix score)*2 \n"
      << "        --gap-extend or -e  gap extend cost, default min(matrix score)\n"
      << "        -o output file\n"
      << "        --self seqfile of multiple sequences. Do pairwise compare.\n";
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
void alignSimple(const bioseq &s1, const bioseq &s2, const int seq1begin, 
      const int seq2begin, ostream &ous, const Progparam &param);
int alignSimpleSelf(const string &seqfile, const Progparam &param);

int main(int argc, char *argv[]) {
   int i = 1;
   string file1, file2, outfile;
   int reverseComplement = 0;
   int seq1begin=1, seq1end=-1, seq2begin=1, seq2end=-1;
   int gapOpen = 2;
   int gapExtend = 1;
   Progparam param;
   bool simpleMethod = true;
   bool SELF_ALIGN=false;
   while (i < argc) {
      if (string(argv[i]) == "-1") { file1=argv[++i]; }
      else if (string(argv[i]) == "-2") { file2=argv[++i]; }
      else if (string(argv[i]) == "-o") { outfile=argv[++i]; }
      else if (string(argv[i]) == "-r") { reverseComplement=atoi(argv[++i]); }
      else if (string(argv[i]) == "-b1") { seq1begin=atoi(argv[++i]); }
      else if (string(argv[i]) == "-e1") { seq1end=atoi(argv[++i]); }
      else if (string(argv[i]) == "-b2") { seq2begin=atoi(argv[++i]); }
      else if (string(argv[i]) == "-e2") { seq2end=atoi(argv[++i]); }
      else if (string(argv[i]) == "--gap-open" 
            || string(argv[i]) == "-g") { 
         gapOpen=atoi(argv[++i]); 
         param.g=gapOpen;
      }
      else if (string(argv[i]) == "--gap-extend"
            || string(argv[i]) == "-e") { 
         gapExtend=atoi(argv[++i]); 
         param.e=gapExtend;
      }
      else if (string(argv[i]) == "--match-score"
            || string(argv[i]) == "-m") { param.m=atoi(argv[++i]); }
      else if (string(argv[i]) == "--mismatch-score"
            || string(argv[i]) == "-s") { param.s=atoi(argv[++i]); }
      else if (string(argv[i]) == "--self") { 
         SELF_ALIGN=true;
         file1=string(argv[++i]); 
      }
      else {
         file1 = argv[i];
         if (i+1 < argc && argv[i+1][0] != '-') {
            file2 = argv[++i];
         }
      }
      ++i;
   }

   if (SELF_ALIGN) {
      return alignSimpleSelf(file1, param);
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
      outfile += ".global.aln";
   }
   ofstream ouf(outfile.c_str());
   if (!ouf) {
      cerr << "Failed to open " << outfile << " for alignment output!\n";
      return 1;
   }

   bioseq seq1, seq2;
   seq1.read(file1);
   seq2.read(file2);
   shrinkSequence(seq1, seq1begin, seq1end);
   shrinkSequence(seq2, seq2begin, seq2end);

   //cerr << "just read 2 sequences\n";
   //cout << seq1.getName() << endl << seq2.getName() << endl;
   //Matrix scoreMatrix;
   //scoreMatrix.setPath("home/zhouke/src/proj/seqaln/matrix");
   // align DNA
   if (seq1.guessType() == DNASEQ) {
      DNA dna1(seq1);
      DNA dna2(seq2);
      if (reverseComplement == 1) dna1.revcomp();
      else if (reverseComplement == 2) dna2.revcomp();
      else if (reverseComplement > 0) {
         cerr << "you can only use 1 or 2 for the -r option\n";
      }
      if (simpleMethod) 
         alignSimple(dna1, dna2, seq1begin, seq2begin, ouf, param);
      else 
         alignDNA(gapOpen, gapExtend, dna1, dna2, seq1begin, seq2begin, ouf);
   }
   else { // protein align
      if (simpleMethod) 
         alignSimple(seq1, seq2, seq1begin, seq2begin, ouf, param);
      else 
         alignProtein(gapOpen, gapExtend, seq1, seq2, seq1begin, seq2begin, ouf);
   }
   cout << "aligment written to file: " << outfile << endl;

   return 0;
}

void alignSimple(const bioseq &s1, const bioseq &s2,
      const int seq1begin, const int seq2begin, ostream &ous,
      const Progparam &param) {
   SimpleScoreMethod sm(param.m, param.s, param.g, param.e);
   Dynaln<SimpleScoreMethod> aln(s1, s2);
   aln.setMatrix(sm);
   aln.runglobal(seq1begin-1, seq2begin-1);
   aln.printAlign(ous);
   showAlignInfo<SimpleScoreMethod>(aln, cout);
}

int alignSimpleSelf(const string &seqfile, const Progparam &param) {
   vector<DNA> allseq;
   ifstream inf(seqfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << seqfile << endl;
      exit(1);
   }
   DNA dna;
   while (dna.read(inf)) {
      allseq.push_back(dna);
   }
   cout << allseq.size() << " sequences stored in memory\n";
   string outfile=seqfile;
   if (outfile.rfind('.') != string::npos) {
      outfile = outfile.substr(0,outfile.rfind('.')) + "_self.aln";
   }
   else outfile += "_self.aln";
   // set up output
   ofstream ous(outfile.c_str());
   SimpleScoreMethod sm(param.m, param.s, param.g, param.e);
   Dynaln<SimpleScoreMethod> aln, rcaln;
   aln.setMatrix(sm);
   rcaln.setMatrix(sm);
   float rcCutoff=0.65;
   int backwardCnt=0;
   time_t timer1, timer2;
   time(&timer1);
   for (size_t i=0; i<allseq.size()-1; ++i) {
      for (size_t j=i+1; j<allseq.size(); ++j) {
         //cerr << "working on " << i << "," << j << " ";
         aln.setseq(allseq[i], allseq[j]);
         aln.runglobal();
         if (aln.getIdentity()>0.85 && (aln.getCov1() > 0.5 || aln.getCov2()>0.5)) {
            aln.printAlign(ous);
         }
         else {
            if (aln.getIdentity() < rcCutoff) {
               ++backwardCnt;
               ous << string(20, '-') << " reverse complement " << string(20, '-') << endl;
               dna=allseq[j].revcompCopy();
               rcaln.setseq(allseq[i], dna);
               rcaln.runglobal();
               if (rcaln.getIdentity()>0.85 && (rcaln.getCov1()>0.5 || rcaln.getCov2()>0.5)) {
                  rcaln.printAlign(ous);
               }
               /*
               rcaln.printAlign(ous);
               if (rcaln.getIdentity() > aln.getIdentity()) {
                  allseq[j]=dna;
               }
               */
            }
         }
      }
      //cout << endl;
   }
   time(&timer2);
   int forwardCnt=allseq.size()*(allseq.size()-1)/2;
   double seconds=difftime(timer2,timer1);
   cout << "took " << seconds << " seconds to do "
      << forwardCnt << " forward and "
      << backwardCnt << " backward alignments\n"
      << seconds/(forwardCnt+backwardCnt) << " seconds per align\n";
   cout << "result written to " << outfile << endl;
   return 0;
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
   alnp.runglobal(seq1begin-1, seq2begin-1); // seqXbegin for labeling the alignment
   alnp.printAlign(ous);
   showAlignInfo<ProteinScoreMethod>(alnp, cout);
}

void alignDNA(const int gapo, const int gape, const DNA &dna1, const DNA &dna2,
      const int seq1begin, const int seq2begin, ostream &ous) {
   Dynaln<NucleicScoreMethod> aln(dna1, dna2);
   if (gapo < 1 || gape < 1) {
      aln.setGapParameter(gapo, gape);
   }
   aln.runglobal(seq1begin-1, seq2begin-1);
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
