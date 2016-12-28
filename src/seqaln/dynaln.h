#ifndef DYNALN_H
#define DYNALN_H

#include "bioseq.h"
#include "matrix.h"
#include <list>
#include <vector>
#include <exception>
#include <stdexcept>

namespace orpara {
// this is deprecated and replaced with the template version
// all binaries depends on this one needs to be rewritten
class AlnInputException : public exception {
   private: 
      string message;
   public:
      AlnInputException(const string &msg) throw() : message(msg) { }
      const char* what() const throw() { return message.c_str(); }
      ~AlnInputException() throw() { }
};

/* experimental, might not be needed 
struct Alignedcol {
   Alignedcol() { }
   Alignedcol(const char t, const char m, const char b) 
      : top(t), middle(m), bottom(b) { }
   char top;
   char middle;
   char bottom;
};

struct Alnstring {
   string top;
   string middle;
   string bottom;
};

struct Alnsummary {
   int score;
   int idencnt;
   int simcnt;  // excluding identical residues
   int numgaps1;
   int numgaps2;
   int gaplen1;
   int gaplen2;
};
*/

/** Define the Alignment Type.
 * Two types: GLOBAL and LOCAL
 */
enum ALNTYPE {GLOBAL, LOCAL};
//enum POINTER {diag, top, left}; // for traceback

/** Class Dynaln: Dynamic Alignment Algorithm.
 * sequence input can only be inputed as a pointer to 
 * a bioseq object.  We might need to remove the const
 * constraint, this way we can change the underlying 
 * sequence easily and cheaply.
 *
 * Matrix is a pointer that this object should not
 * be responsible to deallocate.
 *
 * I am using 0-based index. You simply needs to add 1 to 
 * make it 1-based index.
 */
class Dynaln {
   public:
      /**
       * Default constructor.
       * Default gap insert -8, gap extend -2
       */
      Dynaln() : gapi(-8), gape(-2), ST(), seq1(0), seq2(0),
         S(0), M(0), IX(0), IY(0), alnidx(), 
         Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
         idencnt(0), simcnt(0), 
         numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
         topaln(), middle(), bottomaln(), topruler(), bottomruler() { }
      /**
       * Make a dynamic alignment object Dynaln with two input
       * sequences s1 and s2. The sequence type of s1 and s2 
       * could be both protein or nucleic acids, but not mixed;
       * although, I am planning to implement mixed alignments
       * in the future.
       *
       * @param s1 first sequence
       * @param s2 second sequence. The two sequence can have different
       *        cases, for comparision, lower case was converted to upper case.
       *
       * The scoring table ST will use the default parameters 
       * that are determined by the maxtri class. Right now 
       * this is the blosum50 matrix stored at 
       *    /home/kzhou/proj/seqaln/matrix
       *
       * For repeated usage of this object, it is better to avoid
       * constructing this object repeatedly because reading this
       * matrix file is a disk operation which is slow. Please
       * consider using setSeq() method.
       */
      Dynaln(const bioseq &s1, const bioseq &s2) 
         : gapi(-8), gape(-2), ST(), seq1(&s1), seq2(&s2), 
            S(0), M(0), IX(0), IY(0), alnidx(), 
            Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
            idencnt(0), simcnt(0), 
            numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
            topaln(), middle(), bottomaln(), topruler(), bottomruler()
      { 
         // if nucleic acids, then should switch from default blosum50
         // to nucleic acid matrix
         if ( (s1.getSequenceType() & NUCLEIC_ACID & s2.getSequenceType()) > 0)
         ST.setMatrix("NUC.4.4");
      }
      /**
       * Constructor to build a Dynamic alignment object.
       *
       * @param s1 first sequence
       * @param s2 second sequence
       * @param ma score matrix name such as blosum62
       */
      Dynaln(const bioseq &s1, const bioseq &s2, const string &ma) 
         : ST(ma), seq1(&s1), seq2(&s2), 
            S(0), M(0), IX(0), IY(0), alnidx(), 
            Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
            idencnt(0), simcnt(0), 
            numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
            topaln(), middle(), bottomaln(), topruler(), bottomruler()
      { 
         gapi=2*ST.getMinScore(); 
         gape=ST.getMinScore()/3; 
         // forgot why gape [-3, -1]
         if (gape > -1) gape = -1; 
         if (gape < -3) gape = -3; 
      }
      /** It does the alignment, return the index of the 
       * final cell (in terms of the sequence)
       * in the alignment, this can be used for
       * traceback.  The max score will be set.
       *
       * @return the max score of the global alignment.
       */
      int global() throw(AlnInputException);
      /**
       * This method is used if you want the alignment
       * The alignment result will be stored internally.
       */
      int runglobal(const int delta1=0, const int delta2=0) {
         int s=global(); buildResult(delta1, delta2); return s; }
      /** 
       * Run the local algorithm
       */
      int local();

      /** 
       * Run the local algorithm, then build alignment result for output.
       * use buildResult() to build the alignment result.
       */
      int runlocal(const int delta1=0, const int delta2=0) {
         int s=local(); buildResult(delta1, delta2); return s; }
      //void globalTraceBack();
      ~Dynaln();

      /** 
       * This methods build addional information from the
       * initial phase of the algorithm.
       *    1. call traceback.
       *    2. clear existing result (call clearResult())
       *    3. then use the alnidx list to update other information:
       *       idencnt: counts for identical resudues
       *       simcnt: similar residues
       *       numgaps1, numgaps2: number of gaps in seq1 and seq2
       *       gaplen1, gaplen2: total gaps in both sequences.
       *       topaln: aligned string of seq1
       *       bottomaln: aligned string of seq2
       *       middle: alignment middle line (| for identical : for similar)
       * This method create additional redundant information for
       * easy usage and user interface.
       *
       * For global alignment, initial and terminal gaps are not counted.
       * So the alignment begin index may point to the first aligned residue
       * such as
       * 1         11        21 
       * +         +         +         +  
       * CGCTTTGAATGTTAAGAAATCTC
       *                    ||||
       * -------------------TCTC
       * +         +         +         +  
       *                     2  
       *                    topBeginIndex() is 19 not 0
       */
      void buildResult(const int delta1=0, const int delta2=0);
      /** This method will build the alignment information
       * topaln, buttaln, and middle string will be filled
       *
       * This method more generic, derived class don't have
       * to implement this. The linear space version use this 
       * method.
       * @param delta1 shift index for the first sequnce.
       * @param delta2 shift index for the second sequence.
       * The above parameters are used for aligning subsequences
       * where you need to mark the subsequnce even you have used
       * only the subsequence for the initial alignment.
       */
      void buildAlnInfo(const int delta1=0, const int delta2=0);

      ///////// Output Functions //////////////////////////
      /** 
       * output the alignment result in text format
       * To use this function you must call buildResult()
       */
      void printAlign(ostream &ous, const int w=70) const;
      /** 
       * return the alignment top line as string
       */
      const string& getTopAln() const { return topaln; }
      /** return the bottom line of the alignment as text
       */
      const string& getBottomAln() const { return bottomaln; }
      const bioseq& getTopSequence() const { return *seq1; }
      const bioseq& getBottomSequence() const { return *seq2; }
      /**
       * return the consensus sequence, and the number of residues
       * that are different.
       */
      pair<string, int> getNucleicConsensus() const;

      /**
       * @return at most 2 sequences. In case there is 
       * very low error rates, but high indel rates.
       * Only output at most 2 sequences one of them is more likely
       * to be correct.
       */
      pair<string, string> getNucleicConsensus2() const;

      /** return the match as an array of integers 
       * 0: no match, 1: similar, 2: match
       * @param marr the vector container to hold the result.
       *        This is passed by reference. This function will
       *        first empty this vector, then filling it up.
       * @param countsim a flag to tell whether to count 
       *        similar matches or not.
       * @param append will append to the ngmarr if set true.
       *        will not clear the input vector.
       * To get identity you need to divide the final value by 2
       * This function is mainly designed for bootstrap.
       * sum of the numbers/length/2 is the non-gapped identity.
       */
      void getNgMatchArray(vector<int> &ngmarr, bool append=false, bool countsim=true) const;
      /** similar to getNgMatchArray except it counts the gaps
       * -1: gap, 0: no match, 1: similar, 2: match
       *  For computing identity, gap should be counted 0.
       *  We distingush between gap and no match.
       */
      void getMatchArray(vector<int> &marr, bool countsim=true) const;
      /** output the result in summary table format
       * the delimiter is usually TAB or COMMA.
       * The fields are:
       * seq1name, seq1 length, seq2name, seq2length,
       * score, identity count, similar count,
       * alignment lentth, numgaps of seq1, number gap seq2,
       * gaplen seq1, gaplen seq2, seq1begin, seq1end,
       * seq2 begin, seq2end, entropy seq1, entropy seq2
       *
       * The entropy is the single residue entropy of the
       * subsequence involved in the alignment.
       *
       * @param dl  Delimitor, default comma.
       * @param ibase index base default 0. Use 1 for 1-based index
       * @return the output as table format in string
       */
      string toDelimitedString(const string &dl=",", int ibase=0) const;

      /** print out the two input sequence in string format
       * into STDOUT
       * A debug function.
       */
      void showseq() const {
         cout << seq1->toString() << "\n"
            << seq2->toString() << endl; }

      ////////////  Informational Functions ///////////////////////////////
      /**
       * return the score of the alignment
       */
      int getScore() { return Smax; }
      /**
       * @return the alignment length.
       */
      int getAlnlen() const { return alnidx.size(); }
      float getCov1() const { return alnidx.size()/(float)seq1->length(); }
      float getCov2() const { return alnidx.size()/(float)seq2->length(); }
      /**
       * @return the length of sequence1. If not set then throw exception.
       */
      int getSeq1Length() const throw (std::runtime_error);
      /**
       * @return the length of sequence 2.
       * @throws exception of sequence2 is not set or out of scope
       */
      int getSeq2Length() const throw (std::runtime_error);
      /**
       * @return the name of the first sequence.
       */
      string getSeq1Name() const {
         return seq1->getName();
      }
      /**
       * @return the name of the second sequcne.
       */
      string getSeq2Name() const {
         return seq2->getName();
      }
      string getSequence1AsString() const { return seq1->toString(); }
      string getSequence2AsString() const { return seq2->toString(); }
      /**
       * Get the sequence identity as a fraction range from 0 to 1
       */
      float getIdentity() const { return (float)idencnt/getAlnlen(); }
      float getSimilarity() const { return (float)(idencnt+simcnt)/getAlnlen(); }
      /**
       * Exclding the terminal gap when count the identity.
       * Identical/(alignlen - terminal gaps)
       * This is meaningful for global alignment.
       */
      float getNoTerminalGapIdentity() const;
      /** count the number of gaps in the two sequences
       *
       * @return pair of the number of gaps for both sequences
       */
      pair<int,int> numgaps() const { 
         return make_pair(numgaps1, numgaps2); }
      /**
       * @return number of gaps in the first sequence.
       */
      int getNumgaps1() const { return numgaps1; }
      /**
       * @return number of gaps in the second sequence.
       */
      int getNumgaps2() const { return numgaps2; }
      /**
       * @return total length of gap in the first sequence.
       */
      int getGaplen1() const { return gaplen1; }
      /**
       * @return total gap length in the second sequence.
       */
      int getGaplen2() const { return gaplen2; }

      /**
       * The aligned portion has no gap and 100% identity
       */
      bool isPerfectAlign() const { return numgaps1==0 && numgaps2==0 && getAlnlen() == idencnt; }
      /**
       * All of sequence 1 matched perfectly
       */
      bool isPerfectSeq1Align() const {
         return numgaps1==0 && numgaps2==0 && idencnt == seq1->length(); }
      /**
       * All of sequence 2 matched perfectly.
       */
      bool isPerfectSeq2Align() const {
         return numgaps1==0 && numgaps2==0 && idencnt == seq2->length(); }
      //const Matrix& getScoringMatrix() const { return ST; }
      const Matrix* getScoringMatrix() const { return &ST; }
      /** 
       * 0-based start index of the sequence1 in alignment.
       * For global alignment, the terminal gap are not counted.
       * For example:
       * 1         11        21 
       * +         +         +         +
       * CGCTTTGAATGTTAAGAAATCTC
       *                    ||||
       * -------------------TCTC
       * +         +         +         +
       *                     2  
       * has 19 as topBeginIndex()
       */
      int topBeginIndex() const { return seq1begin; }
      /**
       * 0-based start index of sequenc2 in alignment.
       */
      int bottomBeginIndex() const { return seq2begin; }
      /**
       * 0-based end index of sequence1 in alignment.
       * End index in the alignment. Not counting unaligned region in the end
       * region.
       */
      int topEndIndex() const { return seq1end; }
      /**
       * 0-based end index of sequence2 in alignment.
       */
      int bottomEndIndex() const { return seq2end; }

      /**
       * return the index of the first residue in the alignment
       */
      int topAlignBeginIndex() const;
      int bottomAlignBeginIndex() const;

      ///////// Input and Parameter Settings /////////////////////////
      void setMatrix(const string &mat) { 
         ST.setMatrix(mat); 
         setGapParameterFromMatrix();
      }
      /**
       * make a copy of mat
       * @param mat input matrix to be used by the aligner.
       */
      void setMatrix(const Matrix &mat) { 
         ST=mat; setGapParameterFromMatrix();
      }

      /** 
       * memory allocation should be handled here
       * for repetitive calling of the algorithm.
       * Set the first sequence to the given sequence.
       */
      void setSeq1(const bioseq &sq) { seq1 = &sq;  }
      /**
       * Set the second sequence. Internally, we only keep the pointer,
       * the external object should not be outof scope.
       */
      void setSeq2(const bioseq &sq) { seq2 = &sq; }
      /**
       * set both sequences. same as setseq() method,
       * just different spelling.
       */
      void setSeq(const bioseq &sq1, const bioseq &sq2) {
         seq1=&sq1; seq2=&sq2; }
      /** set up the input sequences for comparison.
       * Not Making a copy!
       * For efficiency, the pointer to sq1 and sq2 are stored,
       * so the bioseq object must not be temporary!
       * @param sq1 the first sequence to use for comparison
       * @param sq2 the second sequence for comparision
       *
       * Note: both sequence must not be temporary, this
       *       function only takes the pointer to 
       *       the sequence objects.
       */
      void setseq(const bioseq &sq1, const bioseq &sq2) {
         seq1=&sq1; seq2=&sq2; }
      //void setSeq(const string &sq1, const string &sq2) {
       //  *seq1=sq1; *seq2=sq2; }

      /**
       * Set the gap insert parameter, default was -8
       * Should be an negative number. If given a positive number, I will flip
       * it to negative.
       * The smaller (more negative) the fewer gaps in the alignment.
       */
      void setGapInsert(int ins) { gapi = ins <= 0? ins : -ins; }
      /**
       * @param ext a negative number -2 is the default.
       *   The smaller the shorter of the gaps. The gap extention
       *   panelty is also relative to the scoring matrix.
       *   If the matrix have large scores, then the extention
       *   will be less penalized.
       *
       * In reality, you want large negative numbers for gap insert parameter,
       * and small negative numbers for gap extend.
       * This method will flip a positive number to negative.
       */
      void setGapExtend(int ext) { gape = ext <= 0? ext : -ext; }
      /**
       * set two gap parameters at the same time.
       * @param go gapOpen or gap insert score. Should be a negative number.
       *        if positive numbers given it will be converted to negative 
       *        value.
       * @param ge gap extend cost. A negative value.
       */
      void setGapParameter(const int go, const int ge) {
         setGapInsert(go); setGapExtend(ge); 
      }

      /**
       * Assume that maxtri has been set.
       * Use the gap parameters stored in matrix ST.
       * If the matrix gap parameters are the defaults (0), then
       * this function call will be NO-OP. This function will also 
       * make sure the values are negative because the methods:
       *   setGapExtend() and setGapInsert() both check for negative 
       *   values.
       */
      void setGapParameterFromMatrix();

      void setAllParameter(const Matrix& mat, const int go, const int ge)
      { ST=mat; setGapInsert(go); setGapExtend(ge); }


      /**
       * Headers for tabular output
       * @param dl delimiter string.
       */
      static string headers(const string &dl);

      /** default character for gap in aligned sequences */
      static const char gapchar='-';
      static const char idenchar='|';
      static const char simchar=':';
      static const int simcut=0;
      /**
       * for marking the sequence alignment, every 10 residues
       */
      static const int markevery=10;
      /** use & | operator to figure out which direction
       * store the traceback pointers in an integer.
       * Multiple pointers can be stored as sum of the 
       * individual pointers
       * DIAG+TOP=3 etc...
       */
      static const int PTRNULL=0;
      static const int PTRDIAG=1;
      static const int PTRTOP=2; // Ix
      static const int PTRLEFT=4; // Iy

      ////////// debug functions /////////
      /** 
       * Debug fucntion.
       * This is used after the algorithm is done. Shows the alignment result
       * matrix, not the scoring matrix.
       * */
      void debug_showmatrix(ostream &ous) const;
      /**
       * Show gap parameters and scoring matrix
       */
      void showParameters() const;

   protected:
      /** set gap values */
      void countnumgaps();
      /** this method is used to clear secondary results that
       * are derived from the first phase of alignment (usually
       * method call of global() or local()) so that
       * this object can be used to do repeated alignments
       * from different inputs. For efficiency.
       * 
       * It is called right after the global() or local()
       * call so that the user will not accedentally use
       * previous rounds of results for this round of alignment
       * if multiple rounds of alignment has been going on within
       * this object.
       */
      void clearResult();
      //void allocmemLS();
      
      /* Nr is the number of row (seq1 length + 1)
       * Nc is the number of columns (seq2 length + 1)
       */
      //void allocmem(const int Nr, const int Nc);
      /** allocate needed memory according to the length
       * of the input sequences
       * 
       * For repeated use of this object, if previous
       * allocated memory is large enough then this
       * object will reuse the old memory without reallocate
       * new memory.
       */
      void allocmem();
      //void traceback(int &i, int &j, int *&Sc, bool local=false);
      
      /*
       * the traceback function produce the alignment result
       * and put it into the alnidx list structure.
       * It also set the begin,end of each sequence that.
       * For global alignment, it will traceback the begining
       * gaps. The ending gap will also be removed from the
       * end of the sequence.
       *
       * i, and j are the matrix location index 0-based
       */
      //void traceback(int &i, int &j, bool local=false);
      //void traceback(int &i, int &j);
      /** this version works with the pointers
       * and only affect the alnidx list.
       */
      void tracepointer(int &i, int &j);
      /**
       */
      void findAlnBoundary();
      //void initialize();
      int gapi, gape;  // should be negative numbers
      //Matrix ST; // scoring table, use protein BLOSUM50 as default.
      //We should use pointers
      /**
       * For multiple usage we should use pointers
       * ST also contain alpha (gapi) and beta (gape).
       * For efficiency, I am using this class's value to avoid
       * an extra level of indirection.
       */
      Matrix ST; // scoring table, use protein BLOSUM50 as default.
      // for DNA use NUC.4.4
      /** Only storing the pointer to input sequences! */
      const bioseq *seq1;
      const bioseq *seq2;
      // the code were assigned when allocating memory

      /**
       * C1 is sequece1->getcode() 
       *    C1=seq1->getcode();
       *    getcode is a virtual function, so bioseq and DNA
       *    use different algorithms.
       * needs to be updated properly when sequence got updated
       */
      const int* C1; 
      /**
       * C2 is sequence2->getcode()
       * Same as C1: code for the sequence.
       */
      const int* C2;
      /** simulated 2-D array of integers
       *  size (seq1.length()+1)*(seq2.length()+1)
       *  For internal use only.
       *  Will store the trace-back pointer.
       */
      int* S; 
      /** the size of the S 
       *  for the next round, memorize the last state
       *  Will not down-size the matrix S for repeated use.
       * */
      int Ssize;
      int *M, *IX, *IY; // array for computing the scores
      //int* Itop;  // for temporary usage
      int numcol;  // for the next round
      //int* P; // for trace back pointer
      //int Ileft;  // for temporary usage
      //stack<pair<int, int> > alnidx;
      
      /** alnidx contains final alignment result.
       * This is the primary result of the alignment algorithm.
       *
       * the numbers are the index in the two sequences.
       * -1 for a gap
       * we can print this result in many different ways
       */
      list<pair<int, int> > alnidx;
      //int score;  // assigned in the trace back step.
      int Smax;  // assigned before the trace back step.
      // in local it is use during the build up phase
      // in global it is the last cell m,n
      /** the two numbers are use to record the
       * traceback starting point in terms of index
       * in the input sequences.
       * For global, it is seq1len-1, seq2len-1
       */
      int Smaxi, Smaxj;  // for local dynamic square memory
      /** start and end of alignment, if local then
       * it is the aligned part, if global, it 
       * excludes the begining and ending gaps
       * 0-based index, inclusive [b,e]
      */
      int seq1begin, seq1end, seq2begin, seq2end;

      // alignment summary and detailed information
      // The following members are for holding secondary results
      /**
       * After computing the alignment results the following will be avaialbe
       * for query
       * Identical residues
       */
      int idencnt;
      /**
       * Similar residues (including identicals)
       */
      int simcnt;
      int numgaps1;
      int numgaps2;
      int gaplen1;
      int gaplen2;

      //int alnlen;
      //vector<Alignedcol> alnresult;
      /** string representation of the alignment
       * topaln    ABCDEFG
       * middle    |: |:
       * bottomaln AEXDDWM
       * The middle line use space for both mismatch and gap
       * so it is not very useful for other purpose other than
       * printing the sequence alignment in plain text format.
       */
      string topaln, middle, bottomaln;  // aligned sequences
      string topruler, bottomruler;
      ALNTYPE alntype;
};

/** Class LSDynaln
 *
 * Linear space Dynamic Alignment.  This is a derived class
 * of the Dynaln. The memory allocation is different so that
 * long sequences can be aligned. 
 *
 * Uses a kind of divided and conquer method.
 *
 */
class LSDynaln : public Dynaln {
   public:
      LSDynaln() : Dynaln() { }
      LSDynaln(const bioseq &s1, const bioseq &s2) 
         : Dynaln(s1,s2) { }
      LSDynaln(const bioseq &s1, const bioseq &s2, const string &mat) 
         : Dynaln(s1,s2,mat) { }
      ~LSDynaln();
      /** the linear space version
       * b1: begin index of sequence1 (0-based index)
       * e1: one-passed the end index of sequence1 (0-based index)
       * @return the best score
       */
      int global() throw(AlnInputException);
      int local();
      /** 
       * Overwrite parent class method. Specific to linear space
       */
      int runlocal(const int delta1=0, const int delta2=0) {
         int s=local(); buildResult(delta1, delta2); return s; 
      }
      //void global();
      /* globalLS in reverse direction
       */
      //int globalR(int b1, int e1, int b2, int e2);
      
      /* b1, e1, b2, and e2 are 0-based index of the 
       * first and the second sequences respectively.
       * if (b1<e1 and b2<e2)
       * then it does forward computation and store
       * the result in S.
       * when b1>e1 and b2>e2 then it does reverse
       * computation and store the result in SR
       *
       * return a pair of pointers to the top and bottom
       *    arrays respectively.
       *    This version is incorrect in detail.
       */
      //pair< int*, int*> computeScore(int b1, int e1, int b2, int e2);
      /* return the start index of the top and bottom rows
       * in the four matrices
       * This can be used to retrieve the top/bottom rows of the
       * three scoring matrices: M, IX, IY
       */
      pair<int, int> computeScoreFW(int b1, int e1, int b2, int e2);
      /** works on the reverse matrices
       * b1>e1 | b2>e2
       */
      pair<int, int> computeScoreBW(int b1, int e1, int b2, int e2);

      /** arguments:
         * 0-based index of sequences
         *  b1, e1:  begin and end of sequence 1
         *  b2, e2: begin and end of sequence 2
         *  Here b1<e1 and b2<e2
         *
         *  The result will be consolidated in the alnidx list.
       */
      int path(int b1, int e1, int b2, int e2);

      /* Completly different from the full dynamic program traceback.
       * Only the last two rows of the result matrix are in memory!
       *
       * given the two arrays M.first(TOP row), M.second (BUTTOM row)
       * and the 0-based index r, c in terms of the matrix
       * len1 is the length of the range of sequence 1
       * len2 is that of the sequence 2
       * if len1 and len2 are given, then it will use the 
       * reverse of the sequence for tracing.
       * Arguments: 
       *    r row index
       *    c column index, r and c specify the starting point
       *    for traceback and it is the cell with the best score.
       */
      //list<pair<int,int> > traceback(pair<int, int> mrow, int b1, int &e1, int b2, int &e2, bool rev=false);
      /** arguments: 
       *    mrow  index of the last two rows of the matrix (top,bottom)
       *    This function will decrement e1 and e2.  If they have become
       *    smaller than b1 and b2 respectively, then the job has been
       *    finished and there is no need to do another recursion call.
       *    The caller (path()) function uses this information.
       * return: path in forward direction
       */
      list<pair<int,int> > tracebackFW(pair<int,int> mrow, int b1, int &e1, int b2, int &e2);
      /**
       * arguments:
       *    mrow is the matrix row, a pair of integers for the top and bottom
       *       array starting point.  This information was used to locate the
       *       starting point for the MR, IXR, and IYR arrays.  The first
       *       integer starts the top array, and the second integer starts the
       *       bottom array.  MR + mrow.first --> top array.
       * Return the path in the forward direction, there is no need to reverse it
       */
      list<pair<int,int> > tracebackBW(pair<int,int> mrow, int b1, int &e1, int b2, int &e2);
      //void showaln();  // testing function

      /** Overwrite the base class version.
       * There are some commonality, need to be factored out
       * in the future.
       */
      void buildResult(const int delta1=0, const int delta2=0);

      void debug_showmatrixForward(pair<int,int> &rows, const int col);
      void debug_showmatrixBackward(pair<int,int> &rows, const int col);

   private:
      /** allocate memory use two rows for top problem
       * two rows for bottom problem
       * S and SR have only one row to store the 
       * score of the last row.
       */
      void allocmem();
      int *MR, *IXR, *IYR;
      //int *SR;  // for LOCAL this is used to memorize the traceback pointer
      //int *ItopR; // not used in local
};
}
#endif
