#ifndef DYNALNT_H
#define DYNALNT_H

// (c) 2002 Kemin Zhou at orpara.com
//File: dynalnt.h

#include "bioseq.h"
#include "scorematrix.h"
#include <list>
#include <vector>
#include <exception>
#include <stdexcept>

/// the following is for class definition
#include <iomanip>
#include <sstream>
#include "strformat.h"
#include <typeinfo>

namespace orpara {
//uncomment if you want to debug the library
//#define DEBUG

/**
 * Exception class, not template.
 */
class AlnInputException : public exception {
   private: 
      string message;
   public:
      AlnInputException(const string &msg) throw() : message(msg) { }
      const char* what() const throw() { return message.c_str(); }
      ~AlnInputException() throw() { }
};

/** 
 * Define the Alignment Type.
 * Two types: GLOBAL and LOCAL
 */
enum ALNTYPE {GLOBAL, LOCAL};

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
template<class T>
class Dynaln {
   public:
      /**
       * Default constructor.
       * Gap parameters will be taken from the matrix.
       * If no parameter is given to the constructor, it will use the default
       * matrix of the given type.
       */
      Dynaln() : ST(), seq1(nullptr), seq2(nullptr),
         S(nullptr), Ssize(0), M(nullptr), IX(nullptr), IY(nullptr), numcol(0), alnidx(), 
         Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
         idencnt(0), simcnt(0), 
         numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
         topaln(), middle(), bottomaln(), topruler(), bottomruler() { setGapParameterFromMatrix(); }

      /**
       * Constructor from a given matrix.
       * @param scmethod  scoring method object one of ScoreMethod class
       *        hierarchy.
       */
      Dynaln(const T &scmethod) 
         : ST(scmethod), seq1(nullptr), seq2(nullptr),
         S(nullptr), Ssize(0), M(nullptr), IX(nullptr), IY(nullptr), alnidx(), 
         Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
         idencnt(0), simcnt(0), 
         numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
         topaln(), middle(), bottomaln(), topruler(), bottomruler() 
      { 
         setGapParameterFromMatrix(); 
#ifdef DEBUG
         cerr << __FILE__ << ":" << __func__ << ":INFO using matrix\n";
         scmethod.show(cerr);
#endif
      }

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
         : ST(), seq1(&s1), seq2(&s2), 
            S(nullptr), Ssize(0), M(nullptr), IX(nullptr), IY(nullptr), alnidx(), 
            Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
            idencnt(0), simcnt(0), 
            numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
            topaln(), middle(), bottomaln(), topruler(), bottomruler()
      { 
         setGapParameterFromMatrix();
         if (!validInput()) throw AlnInputException("Matrix and Input sequence type mismatch");
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
            S(nullptr), Ssize(0), M(nullptr), IX(nullptr), IY(nullptr), alnidx(), 
            Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
            idencnt(0), simcnt(0), 
            numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
            topaln(), middle(), bottomaln(), topruler(), bottomruler()
      { 
         setGapParameterFromMatrix();
         if (!validInput()) throw AlnInputException("Matrix and Input sequence type mismatch");
      }
      /**
       * Constructor from all known resources.
       * @param s1 a const reference to first sequence.
       * @param s2 const reference to second sequence.
       * @param sm scoring matrix used by this aligner.
       */
      Dynaln(const bioseq &s1, const bioseq &s2, const T &sm) 
         : ST(sm), seq1(&s1), seq2(&s2), 
            S(nullptr), Ssize(0), M(nullptr), IX(nullptr), IY(nullptr), alnidx(), 
            Smax(0), seq1begin(-1), seq1end(-1), seq2begin(-1), seq2end(-1),
            idencnt(0), simcnt(0), 
            numgaps1(0), numgaps2(0), gaplen1(0), gaplen2(0),
            topaln(), middle(), bottomaln(), topruler(), bottomruler()
      { 
         setGapParameterFromMatrix();
         if (!validInput()) throw AlnInputException("Matrix and Input sequence type mismatch");
      }

      Dynaln(const Dynaln &da)=delete;
      Dynaln& operator=(const Dynaln &da)=delete;
      /**
       * This should be used in very liminted numbers.
       * It will consume a lots of memory!
       */
      Dynaln(Dynaln &&o) : gapi(o.gapi), gape(o.gape),
         ST(std::move(o.ST)), seq1(o.seq1), seq2(o.seq2),
         C1(o.C1), C2(o.C2), S(o.S), Ssize(o.Ssize),
         M(o.M), IX(o.IX), IY(o.IY), numcol(o.numcol),
         alnidx(std::move(o.alnidx)), Smax(o.Smax),
         Smaxi(o.Smaxi), Smaxj(o.Smaxj),
         seq1begin(o.seq1begin), seq1end(o.seq1end),
         seq2begin(o.seq2begin), seq2end(o.seq2end),
         idencnt(o.idencnt), simcnt(o.simcnt),
         numgaps1(o.numgaps1), numgaps2(o.numgaps2),
         gaplen1(o.gaplen1), gaplen2(o.gaplen2),
         topaln(std::move(o.topaln)), middle(std::move(o.middle)), 
         bottomaln(std::move(o.bottomaln)), topruler(std::move(o.topruler)),
         bottomruler(std::move(o.bottomruler)), alntype(o.alntype)
      {
         o.C1=nullptr; o.C2=nullptr; o.S=nullptr;
         o.M=nullptr; o.IX=nullptr; o.IY=nullptr;
      }
      /**
       * Use with care.
       */
      Dynaln& operator=(Dynaln&& o) {
         if (this != &o) {
            gapi=o.gapi; gape=o.gape; ST=std::move(o.ST);
            seq1=o.seq1; seq2=o.seq2;
            if (C1 != nullptr) delete[] C1; 
            if (C2 != nullptr) delete C2; 
            if (S != nullptr) delete S;
            C1=o.C1; C2=o.C2; S=o.S; 
            o.C1=o.C2=o.S=nullptr;
            Ssize=o.Ssize;
            if (M != nullptr) delete M; 
            if (IX != nullptr) delete IX;
            if (IY != nullptr) delete IY;
            M=o.M; IX=o.IX; IY=o.IY; 
            o.M=o.IX=o.IY=0;
            numcol=o.numcol;
            alnidx=std::move(o.alnidx); Smax=o.Smax;
            Smaxi=o.Smaxi; Smaxj=o.Smaxj;
            seq1begin=o.seq1begin; seq1end=o.seq1end;
            seq2begin=o.seq2begin; seq2end=o.seq2end;
            idencnt=o.idencnt; simcnt=o.simcnt;
            numgaps1=o.numgaps1; numgaps2=o.numgaps2;
            gaplen1=o.gaplen1; gaplen2=o.gaplen2;
            topaln=std::move(o.topaln); middle=std::move(o.middle);
            bottomaln=std::move(o.bottomaln); topruler=std::move(o.topruler);
            bottomruler=std::move(o.bottomruler); alntype=o.alntype;
         }
         return *this;
      }

      /** It does the alignment, return the index of the 
       * final cell (in terms of the sequence)
       * in the alignment, this can be used for
       * traceback.  The max score will be set.
       *
       * @return the max score of the global alignment.
       */
      int global();
      /**
       * This method is used if you want the alignment
       * The alignment result will be stored internally.
       */
      int runglobal(const int delta1=0, const int delta2=0) {
         global(); buildResult(delta1, delta2); return this->Smax; }
      /** 
       * Run the local algorithm
       */
      int local();

      /** 
       * Run the local algorithm, then build alignment result for output.
       * use buildResult() to build the alignment result.
       * @param delta1 use this one for labeling the sequence position
       *     in the first input sequence. Starting position of the
       *     input sequence used for alignment.
       * @param delta2 mark position in the second sequence.
       * You could have used a subsequence for the alignment,
       * then you want to label the positions in the original sequence.
       * In this case you add the delta to the labeling.
       *
       * Note when constructing thread object, the two default
       * arguments must be provided.
       */
      int runlocal(const int delta1=0, const int delta2=0) {
               local(); buildResult(delta1, delta2); return this->Smax; 
      }
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
       * @param delta1 the marking of the position of the first sequence.
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
       * the runlocal and runglobal methods will call the
       * buildResult() method for you.
       * @param w the window width. Default 70 char.
       */
      void printAlign(ostream &ous, const int w=70) const;
      string alignToString(const int w=70) const;

      /**
       * output headers.
       */
      static ostream& printSummaryHeader(ostream &ous, const string &dl="\t",
            bool wantEntropy=false);
      /**
       * Print the summary in tabular format.
       * columns:
       * seq1name seq1len seq2name seq2len score idencnt simcnt
       * alnlen numgap1 numgap2 gap1len gap2len 
       * seq1begin seq1end seq2begin seq2 end
       * optional: seq1entropy_of_aligned_region
       *           seq2entropy_of_aligned_region
       * seq1begin, seq1end, ..., used 1-based index in the sequence.
       */
      ostream& printSummary(ostream &ous, const string &dl="\t", bool wantEntropy=false) const;
      /** 
       * @return the alignment top line as string.
       */
      const string& getTopAln() const { return topaln; }
      /**
       * Middle identical or similar or space.
       *  idenchar='|'; simchar=':';
       */
      const string& getMiddleAln() const { return middle; }
      /** 
       * @return the bottom line of the alignment as text which contains
       *    the aligned residue with potential gat char (default - )
       */
      const string& getBottomAln() const { return bottomaln; }
      /**
       * @return a reference to a const bioseq object
       *   represented by the first sequence. 
       */
      const bioseq& getTopSequence() const { return *seq1; }
      const bioseq* getSeq1() const {
         return seq1;
      }
      const bioseq* getSeq2() const {
         return seq2;
      }
      /**
       * @return a reference to a const bioseq object
       *   represented by the second sequence. 
       */
      const bioseq& getBottomSequence() const { return *seq2; }
      /**
       * return the consensus sequence, and the number of residues
       * that are different.
       * This method only applies to nucleic acid sequences.
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
       * @return the score of the alignment
       */
      int getScore() const { return Smax; }
      /**
       * @return the alignment length.
       */
      int getAlnlen() const { return alnidx.size(); }
      /**
       * This is the part of first sequence that are aligned
       */
      float getCov1() const { return getSeq1AlignedLength()/(float)seq1->length(); }
      /**
       * Fraction of second sequence inside alignment
       */
      float getCov2() const { return getSeq2AlignedLength()/(float)seq2->length(); }
      /**
       * @return the length of sequence1. If not set then throw exception.
       */
      int getSeq1Length() const;
      /**
       * @return the length of sequence 2.
       * @throws exception of sequence2 is not set or out of scope
       */
      int getSeq2Length() const;
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
      char getSeq1charAt(const int i) const {
         if (i >= seq1->length()) 
            throw runtime_error("index " + to_string(i) + " bigger than seq1 len");
         return (*seq1)[i];
      }
      char getSeq2charAt(const int i) const {
         if (i >= seq2->length()) 
            throw runtime_error("index " + to_string(i) + " bigger than seq2 len");
         return (*seq2)[i];
      }
      string getSeq1Substring(int b, int len) const {
         return seq1->substring(b,len);
      }
      string getSeq2Substring(int b, int len) const {
         return seq2->substring(b,len);
      }
      /**
       * @return the number of identical residues in the alighment.
       */
      int getIdentical() const { return idencnt; }
      /**
       * Get the sequence identity as a fraction range from 0 to 1
       */
      float getIdentity() const { 
         if (getAlnlen() == 0) return 0;
         return (float)idencnt/getAlnlen(); 
      }
      /**
       * Similar residues are also counted.
       */
      float getSimilarity() const { 
         if (getAlnlen() == 0) return 0;
         return (float)(idencnt+simcnt)/getAlnlen(); 
      }
      /**
       * Exclding the terminal gap when count the identity.
       * Identical/(alignlen - terminal gaps)
       * This is meaningful for global alignment.
       */
      float getNoTerminalGapIdentity() const;
      /**
       * @return identity excluding gaps.
       */
      float getNogapIdentity() const {
         return float(idencnt+numgaps1+numgaps2)/getAlnlen(); }
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
       * @return the total number of gaps for both sequences.
       */
      pair<int,int> gaplen() const { return make_pair(gaplen1, gaplen2); }
      /**
       * @return total length of gap in the first sequence.
       */
      int getGaplen1() const { return gaplen1; }
      /**
       * @return total gap length in the second sequence.
       */
      int getGaplen2() const { return gaplen2; }
      /**
       * Edit distance is the alignment length - identical count
       */
      int getEditDistance() const {
         return getAlnlen() - getIdentical();
      }

      /**
       * The aligned portion has no gap and 100% identity
       */
      bool isPerfectAlign() const { return numgaps1==0 && numgaps2==0 && getAlnlen() == idencnt; }
      /**
       * All of sequence 1 matched perfectly
       */
      bool isPerfectSeq1Align() const {
         return numgaps1==0 && numgaps2==0 && idencnt == int(seq1->length()); }
      /**
       * All of sequence 2 matched perfectly.
       */
      bool isPerfectSeq2Align() const {
         return numgaps1==0 && numgaps2==0 && idencnt == seq2->length(); }
      static string headers(const string &dl);
      /**
       * @return the pair of aligned positions.
       *         -1 means a gap in the sequence.
       *         each pair is <idx1, idx2> 0-based index position
       *         of a matched pair of residues.
       */
      const list<pair<int, int> >& getAlnindex() const { return alnidx; }
      /**
       * @return the alignment detail in Cigar format as defined in Samtools
       *   using seq1 as reference.
       */
      vector<pair<char,int>> getCigar1() const;
      /**
       * @return the alignment detail in Cigar format as defined in Samtools
       *   using seq2 as reference.
       */
      vector<pair<char,int>> getCigar2() const;
      vector<pair<int, int> > getAlnindexVector() const {
         return vector<pair<int,int> >(alnidx.begin(), alnidx.end());
      }
      /** 
       * @return a reference to the underlying scoring matrix
       */
      const T* getScoringMatrix() const { return &ST; }
      /** 
       * 0-based start index in sequence1 for the alignment.
       * For global alignment, the terminal gap are not counted.
       *
       * \verbatim
       * For example:
       * 1         11        21 
       * +         +         +         +
       * CGCTTTGAATGTTAAGAAATCTC
       *                    ||||
       * -------------------TCTC
       * +         +         +         +
       *                     2  
       * has 19 as topBeginIndex()
       *
       * 6         16        26        36        46        56        66 
       * +         +         +         +         +         +         +
       * ATGCACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGGAGTCCCCTTCTTCTCAT
       * |||||||||||||||||| ||||||||||||||||||||||| ||| |||| || ||||||||||| || |||||||| |
       * ATGCACGGTGTTGACTGACTTCAAGACCTGGCTCCAGTCCAAACTCTTGCCACGGTTGCCGGGAGTTCCTTTCTTCTCGT
       * +         +         +         +         +         +         + 
       * 1         11        21        31        41        51        61  
       * \endverbatim
       *
       * The above topBeginIndex() is 5 (6-1; zero based index), the 
       * bottomBeginIndex() is 0;
       * In the printout for humans I used 1-based index.
       * @return index in sequence1 for the first aligned residue.
       */
      int topBeginIndex() const { return seq1begin; }
      /**
       * 0-based start index of sequence2 in alignment.
       */
      int bottomBeginIndex() const { return seq2begin; }
      /**
       * 0-based end index of sequence1 included in the alignment.
       * End index in the alignment. Not counting unaligned region in the end
       * region.
       */
      int topEndIndex() const { return seq1end; }
      /**
       * 0-based end index of sequence2 included in the alignment.
       */
      int bottomEndIndex() const { return seq2end; }
      int getSeq1AlignedLength() const { return topEndIndex() - topBeginIndex() + 1; }
      int getSeq2AlignedLength() const { return bottomEndIndex() - bottomBeginIndex() + 1; }

      /**
       * return the index of the first residue in the alignment
       */
      int topAlignBeginIndex() const;
      int bottomAlignBeginIndex() const;
      /**
       * return the pair of residues according to the positon of the top
       * sequence. 
       * @param tpos, top 1-based index position.
       * @param bpos will be set to the bottom position 1-based index.
       *         The position for the bottom sequence could be a gap char
       *         if the return value.second == gapchar.
       * @return the matched or mismated residue at the same positon in the
       *         alignment.
       */
      pair<char, char> getResiduesByTopPosition(int tpos, int &bpos);
      /**
       * Given the position of the top sequence, it extract
       * the bottom aligned string. Not the aligned position of the
       * top aligned sequence.
       * @param tpos1 is the top sequence position 1-based
       * @param tpos2 is the end position [tpos1, tpos2] closed range.
       * @return bottom alignment string based on the top sequence
       *   position 1-based.
       */
      string getBottomAlnByTopPosition(int tpos1, int tpos2) const;

      ///////// Input and Parameter Settings /////////////////////////
      /**
       * Not sure this should be use frequently.
       * The type of matrix must match.
       * Replace the underlying matrix with a new matrix named 'mat'.
       */
      void setMatrix(const string &mat) { 
         ST.use(mat); 
         setGapParameterFromMatrix();
      }

      /**
       * make a copy of mat
       * @param mat input matrix to be used by the aligner.
       */
      void setMatrix(const T &mat) { 
         ST=mat; setGapParameterFromMatrix();
      }

      /** 
       * memory allocation should be handled here
       * for repetitive calling of the algorithm.
       * Set the first sequence to the given sequence.
       * I am not checking the type match of sq and ST!
       * For repeated checking of large number of sequence this is too
       * expensive.
       */
      void setSeq1(const bioseq &sq) { seq1 = &sq;  }
      /**
       * Set the second sequence. Internally, we only keep the pointer,
       * the external object should not be outof scope.
       */
      void setSeq2(const bioseq &sq) { seq2 = &sq; }
      /**
       * Alias for setseq()
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

      /**
       * The matrix must be the same type.
       */
      void setAllParameter(const T& mat, const int go, const int ge)
      { ST=mat; setGapInsert(go); setGapExtend(ge); }

      //// for verification /////
      /**
       * @return true if the matrix and the input sequence type match.
       *      otherwise return false.
       */
      bool validInput() const {
         if (!strcmp(typeid(ST).name(), "ScoreMethod") || !strcmp(typeid(ST).name(), "SimpleScoreMethod"))
            return true;
         if (!strcmp(typeid(ST).name(), "NucleicScoreMethod") && seq1->getSequenceType() != NUCLEIC_ACID)
            return false;
         if (!strcmp(typeid(ST).name(), "ProteinScoreMethod") && seq1->getSequenceType() != PROTEINSEQ)
            return false;
         return true;
      }

      /////////////////////////////////////////////////////////
      /** default character for gap in aligned sequences */
      //static const char gapchar='-';
      static const char gapchar;
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
      
      /** allocate needed memory according to the length
       * of the input sequences
       * 
       * For repeated use of this object, if previous
       * allocated memory is large enough then this
       * object will reuse the old memory without reallocate
       * new memory.
       */
      void allocmem();

      /** 
       * This version works with the pointers
       * and only affect the alnidx list.
       * The alnidx list will be updated.
       */
      void tracepointer(int &i, int &j);
      /**
       * Update the begin and end index of the 
       * aligned portion of the two input sequences.
       * If there is no alignment, then all four members will be set to 
       * negative 1 
       * seq1begin=-1; seq2begin=-1; seq1end=-1; seq2end=-1;
       */
      void findAlnBoundary();

      /**
       * Scoring look up table.
       * This is now a parameter for this class.
       * Any of the ScoreMethod class will work.
       */
      T ST; 
      /**
       * Gap parameters for this program.
       * For performance, we have one extra copy here.
       * The ST object also carries this value.
       * should be negative numbers
       */
      int gapi, gape;  

      /** 
       * For performance, only storing the pointer to input sequences! 
       * So the user has to be aware of this. This means that the pointer is
       * only meaningful in the bioseq were in scope.
       * Note: For bioseq class, I am still using lower case first letter,
       * should use upper case for user classes. 
       * */
      const bioseq *seq1;
      const bioseq *seq2;

      // the code were assigned when allocating memory
      /**
       * C1 is sequece1->getcode() 
       *    C1=seq1->getcode();
       *    getcode is a virtual function, so bioseq and DNA
       *    use different algorithms.
       * needs to be updated properly when sequence got updated
       * Note: may  cause race condition in multithreaded code
       *   if they share the same sequence. Call get code first
       *   before sharing the object.
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
       *  This will save the time for reallocating memory unless new requests
       *  are larger than previous ones.
       *  The extra row, column is for boundary condition.
       */
      int* S; 
      /** the size of the S 
       *  for the next round, memorize the last state
       *  Will not down-size the matrix S for repeated use.
       * */
      int Ssize;
      /**
       *  Arrays for computing the scores.
       *  These are all mutable.
       *  Their size should be numcol+1.
       */
      int *M, *IX, *IY; 
      /**
       * remember the number of columns should be
       * the same as last input sequence. 
       * For repeated runs.
       */
      int numcol; 
      
      /** 
       * alnidx contains the final alignment result.
       * This is the primary result of the alignment algorithm.
       *
       * The numbers are the indices in the two sequences.
       * -1 for a gap.
       * we can print this result in many different ways
       * Many secondary results are computed from this data structure.
       */
      list<pair<int, int> > alnidx;
      /**
       * variable used by the algorimth to remember the maximum score
       * after the algorithm is finished.
       */
      int Smax;  
      // in local it is use during the build up phase
      // in global it is the last cell m,n
      /** the two numbers are use to record the
       * traceback starting point in terms of index
       * in the input sequences.
       * For global, it is seq1len-1, seq2len-1, but we set it to make uniform
       * algorithm.
       * This applicable only to square memory allocation, not linear.
       */
      int Smaxi, Smaxj;  
      /** start and end of alignment, if local then
       * it is the aligned part, if global, it 
       * excludes the begining and ending gaps
       * 0-based index, inclusive [b,e]
       *
       * If no match is found then these values will all be
       * set to -1.
      */
      int seq1begin, seq1end, seq2begin, seq2end;

      // alignment summary and detailed information
      // The following members are for holding secondary results
      // Recomputing them will be costly, so we remember them.
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
      /**
       * Number of gaps in the alignment for the first sequence.
       */
      int numgaps1;
      int numgaps2;
      int gaplen1;
      int gaplen2;

      /** 
       * string representation of the alignment
       * \verbatim
       * topaln    ABCDEFG-D
       * middle    |: |:   |
       * bottomaln AEXDDWMKD
       * \endverbatim
       * The middle line use space for both mismatch and gap
       * so it is not very useful for other purpose other than
       * printing the sequence alignment in plain text format.
       */
      string topaln, middle, bottomaln;  // aligned sequences
      /**
       * Delta1 and 2 affects the labeling of the ruler
       */
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
template<class T>
class LSDynaln : public Dynaln<T> {
   public:
      LSDynaln() : Dynaln<T>() { }
      LSDynaln(const bioseq &s1, const bioseq &s2) 
         : Dynaln<T>(s1,s2) { }
      LSDynaln(const bioseq &s1, const bioseq &s2, const string &mat) 
         : Dynaln<T>(s1, s2, mat) { }
      ~LSDynaln();
      /** the linear space version
       * b1: begin index of sequence1 (0-based index)
       * e1: one-passed the end index of sequence1 (0-based index)
       * @return the best score
       */
      int global();
      int local();
      /** 
       * Overwrite parent class method. Specific to linear space
       */
      int runlocal(const int delta1=0, const int delta2=0) {
         local(); buildResult(delta1, delta2); return this->Smax; 
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

//////////// class definition ///////////////////////////


//#define DEBUG

// need to define const members
template<class T> const char Dynaln<T>::gapchar='-';

template<class T>
Dynaln<T>::~Dynaln<T>() {
   if (S != 0) delete[] S; 
   delete[] M;
   delete[] IX;
   delete[] IY;
}

template<class T>
void Dynaln<T>::setGapParameterFromMatrix() {
   setGapInsert(ST.getGapInsert());
   setGapExtend(ST.getGapExtend());
}

template<class T>
void Dynaln<T>::allocmem() {
   int Nr=seq1->length();
   int Nc=seq2->length();
   if (S == nullptr) { // allocate from new memory
      Ssize=Nr*Nc;
      numcol=Nc;
      S = new int[Ssize];
      M = new int[Nc+1]; // for match state
      IX = new int[Nc+1]; // for Insertion of X
      IY = new int[Nc+1]; // for insertion of Y
   }
   else {
      if (Ssize < Nr*Nc) { // old memory not enough
         delete[] S;
         Ssize=Nr*Nc;
         S=new int[Ssize];
      }
      if (numcol < Nc) { // column memory not enough
         delete[] M; 
         delete[] IX; 
         delete[] IY;
         M = new int[Nc+1];
         IX = new int[Nc+1];
         IY = new int[Nc+1];
         numcol=Nc;
      }
   }
   // asigne the code of each sequences
   C1=seq1->getcode();
   C2=seq2->getcode();
}

/*
 * use last direction as a guide to break a tie
 * 0---j--->
 * |
 * i
 * |
 *\|/
 */
template<class T>
void Dynaln<T>::tracepointer(int &i, int &j) {
   if (!alnidx.empty()) {
      alnidx.clear();
   }
   // nothing need to be done if no result for local.
   if (Smax == 0 && alntype == LOCAL) return;
   int Nc=seq2->length();
   //int pi, pj;
   int *ptr = S+(i*seq2->length() + j);
   int lastTrace = PTRNULL; // should have only one direction
   while (i>-1 && j>-1) {
      if ( *ptr == PTRDIAG) {
         alnidx.push_front(make_pair(i,j));
         lastTrace = PTRDIAG;
         ptr -= (Nc+1);
         --i; --j;
      }
      else if (*ptr == PTRLEFT) {
         alnidx.push_front(make_pair(-1, j));
         lastTrace = PTRLEFT;
         --ptr;
         --j;
      }
      else if (*ptr == PTRTOP) {
         alnidx.push_front(make_pair(i, -1));
         lastTrace = PTRTOP;
         ptr -= Nc;
         --i;
      }
      else if ( ((*ptr)&PTRDIAG) == PTRDIAG && ((*ptr)&PTRLEFT) == PTRLEFT
            && ((*ptr)&PTRTOP) == PTRTOP) { // all three directions
         if (lastTrace == PTRDIAG) {
            alnidx.push_front(make_pair(i,j));
            ptr -= (Nc+1);
            --i; --j;
         }
         else if (lastTrace == PTRLEFT) {
            alnidx.push_front(make_pair(-1, j));
            --ptr;
            --j;
         }
         else { // PTRTOP
            alnidx.push_front(make_pair(i, -1));
            ptr -= Nc;
            --i;
         }
      }
      else if ( ((*ptr)&PTRDIAG) == PTRDIAG && ((*ptr)&PTRTOP) == PTRTOP) { // both top and diag
         if (lastTrace == PTRDIAG) {
            alnidx.push_front(make_pair(i,j));
            ptr -= (Nc+1);
            --i; --j;
         }
         else if (lastTrace == PTRTOP) {
            alnidx.push_front(make_pair(i, -1));
            ptr -= Nc;
            --i;
         }
         else { // just pick diagnal
            alnidx.push_front(make_pair(i,j));
            lastTrace = PTRDIAG;
            ptr -= (Nc+1);
            --i; --j;
         }
      }
      else if (((*ptr)&PTRDIAG) == PTRDIAG && ((*ptr)&PTRLEFT) == PTRLEFT)  { // both left and diag
         if (lastTrace == PTRDIAG) {
            alnidx.push_front(make_pair(i,j));
            ptr -= (Nc+1);
            --i; --j;
         }
         else if (lastTrace == PTRLEFT) {
            alnidx.push_front(make_pair(-1, j));
            --ptr;
            --j;
         }
         else { // choose diag
            alnidx.push_front(make_pair(i,j));
            lastTrace = PTRDIAG;
            ptr -= (Nc+1);
            --i; --j;
         }
      }
      else if (((*ptr)&PTRLEFT) == PTRLEFT && ((*ptr)&PTRTOP) == PTRTOP) { // both Left and Top
         if (lastTrace == PTRLEFT) {
            alnidx.push_front(make_pair(-1, j));
            --ptr;
            --j;
         }
         else if (lastTrace == PTRTOP) {
            alnidx.push_front(make_pair(i, -1));
            ptr -= Nc;
            --i;
         }
         else { // choose left
            alnidx.push_front(make_pair(-1, j));
            lastTrace = PTRLEFT;
            --ptr;
            --j;
         }
      }
      else { // NULL case
         if ((*ptr) != PTRNULL) {
            cerr << "it should be the null case\n";
         }
         break; // for local alnment
      }
   }
   if (alntype == GLOBAL) {
      if (j == -1 && i>0) {
         while (i>-1) {
            alnidx.push_front(make_pair(i--,-1));
         }
      }
      else if (i==-1 && j>0) {
         while (j>-1) {
            alnidx.push_front(make_pair(-1,j--));
         }
      }
   }
}


template<class T>
void Dynaln<T>::findAlnBoundary() {
   if (alntype == LOCAL) {
      if (S == 0) {
         seq1begin=-1;
         seq2begin=-1;
         seq1end=-1;
         seq2end=-1;
      }
      else {
         seq1begin=alnidx.begin()->first;
         seq2begin=alnidx.begin()->second;
         seq1end=alnidx.rbegin()->first;
         seq2end=alnidx.rbegin()->second;
      }
   }
   else if (alntype == GLOBAL) {
      // For global alignement, the end will be corrected for the terminal gaps.
      seq1begin=0;
      seq2begin=0;
      list<pair<int,int> >::const_iterator it=alnidx.begin();
      while (it != alnidx.end() && it->first == -1) {
         ++seq2begin;
         ++it;
      }
      it=alnidx.begin();
      while (it != alnidx.end() && it->second == -1) {
         ++seq1begin;
         ++it;
      }
      seq1end=seq1->length()-1;
      seq2end=seq2->length()-1;
      list<pair<int,int> >::reverse_iterator li=alnidx.rbegin();
      while (li != alnidx.rend() && li->first == -1) {
         --seq2end;
         ++li;
      }
      li=alnidx.rbegin();
      while (li != alnidx.rend() && li->second == -1) {
         --seq1end;
         ++li;
      }
   }
}

template<class T>
int Dynaln<T>::topAlignBeginIndex() const {
   list<pair<int,int> >::const_iterator it = alnidx.begin();
   while (it != alnidx.end() && it->first == -1) {
      ++it;
   }
   return it->first;
}

template<class T>
int Dynaln<T>::bottomAlignBeginIndex() const {
   list<pair<int,int> >::const_iterator it = alnidx.begin();
   while (it != alnidx.end() && it->second == -1) {
      ++it;
   }
   return it->second;
}

template<class T>
float Dynaln<T>::getNoTerminalGapIdentity() const {
   int e1len = seq1->length() - seq1end - 1;
   int e2len = seq2->length() - seq2end - 1;
   return float(idencnt)/(getAlnlen() - max(seq1begin, seq2begin) - max(e1len, e2len));
}


/* shows the pointer matrix */
template<class T>
void Dynaln<T>::debug_showmatrix(ostream &ous) const {
   int Nc=seq2->length();
   int i,j;
   ous << "Ssize: " << Ssize << endl;
   ous << "\t*";
   for (j=0; j<seq2->length(); j++) ous << (*seq2)[j] << "\t";
   ous << endl;
   for (i=0; i<seq1->length(); i++) {
      ous << (*seq1)[i] << "\t";
      for (j=0; j<seq2->length(); j++) {
         ous << S[i*Nc+j] << "\t";
      }
      ous << endl;
   }
   ous << endl;
}

template<class T>
void Dynaln<T>::showParameters() const {
   cout << "gap insert: " << gapi << " gap extend: " << gape << endl;
   ST.show();
}

// return score
//pair<int,int> Dynaln<T>::global() {
// pointer for traceback (j-i) diag=0, top=1 (Ix), left=-1 (Iy)
// only store the pointers not the scores, too many
template<class T>
int Dynaln<T>::global() {
   if (seq1->length() < 1) {
      throw AlnInputException("seq1 is empty");
   }
   else if (seq2->length() < 1) {
      throw AlnInputException("seq2 is empty");
   }

   alntype=GLOBAL;
   allocmem();
   int i,j;
   int s1len=seq1->length();  // for short typing
   int s2len=seq2->length();

   M[0]=IX[0]=IY[0]=0;
   M[1]=IX[1]=IY[1]=gapi;
   for (j=2; j<s2len+1; j++) {
      IX[j]=IY[j]=M[j]=M[j-1]+gape;
   }
   int leftM=gapi; 
   int leftIx=gapi;
   int leftIy=gapi;
   int CM=gapi, CIx=gapi, CIy=gapi, backptr; // current values
   int maxDiag, maxTop, maxLeft;
   
#ifdef DEBUG
   cout << "Global alignment matrix:\n";
#endif
   for (i=0; i<s1len; i++) {
#ifdef DEBUG
      cout << "row " << i << endl << leftM << "," << leftIx << "," << leftIy << " | ";
#endif
      for (j=0; j<s2len; j++) {
         maxDiag = max(M[j], max(IX[j], IY[j]));
         CM=maxDiag + ST.lookup(C1[i], C2[j]);
         CIx=max(M[j+1]+gapi, IX[j+1]+gape);
         CIy=max(leftM+gapi, leftIy+gape);
#ifdef DEBUG
         cout << ST.lookup(C1[i], C2[j]) << " " << (*seq1)[i] << "x" << (*seq2)[j] << " "
            << CM << "," << CIx << "," << CIy << " | ";
#endif
         // store pointer
         //Smax=max(CM, max(CIx, CIy));
         if (CM > CIx && CM > CIy) {
            Smax = CM;
            backptr = PTRDIAG;
         }
         else if (CIx > CM && CIx > CIy) {
            Smax = CIx;
            backptr = PTRTOP;
         }
         else if (CIy > CIx && CIy > CM) {
            Smax = CIy;
            backptr = PTRLEFT;
         }
         else if (CM == CIx && CM > CIy) {
            Smax = CM;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            if (maxDiag > maxTop) backptr = PTRDIAG;
            else if (maxTop > maxDiag) backptr = PTRTOP;
            else {
               backptr = (PTRTOP | PTRDIAG);
            }
         }
         else if (CM == CIy && CM > CIx) {
            Smax = CIy;
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxDiag > maxLeft) backptr = PTRDIAG;
            else if (maxLeft > maxDiag) backptr = PTRLEFT;
            else {
               backptr = (PTRLEFT|PTRDIAG); // this is very rare
            }
         }
         else if (CIy == CIx && CIy > CM) {
            Smax = CIy;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxTop > maxLeft) backptr = PTRTOP;
            else if (maxLeft > maxTop) backptr = PTRLEFT;
            else {
               backptr = (PTRLEFT | PTRTOP);
            }
         }
         else { //All three values: CM, CIx, CIy are the same, this is rare!
            Smax = CM;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxDiag > maxTop && maxDiag > maxLeft) backptr= PTRDIAG;
            else if (maxTop > maxDiag && maxTop > maxLeft) backptr= PTRTOP;
            else if (maxLeft > maxDiag && maxLeft > maxTop) backptr= PTRLEFT;
            else if (maxDiag == maxTop && maxDiag > maxLeft) {
               backptr= PTRDIAG | PTRTOP;
            }
            else if (maxDiag == maxLeft && maxDiag > maxTop) {
               backptr= PTRDIAG | PTRLEFT;
            }
            else if (maxTop == maxLeft && maxTop > maxDiag) {
               backptr= PTRTOP | PTRLEFT;
            }
            else {
               //cerr << "it is very rare that the three values are identical: "
               //   << maxDiag << ", " << maxTop << ", " << maxLeft << endl;
               backptr= PTRDIAG | PTRTOP | PTRLEFT;
            }
         }
         //cout << backptr << "\t";
         S[i*s2len+j] = backptr;
         // ready for the next round
         M[j]=leftM; IX[j]=leftIx; IY[j]=leftIy;
         leftM=CM; leftIx=CIx; leftIy=CIy;
      }
#ifdef DEBUG
      cout << endl;
#endif
      M[j]=CM; IX[j]=CIx; IY[j]=CIy;
      leftM=leftIx=leftIy=(M[0] + gape);
   }
   Smaxi=seq1->length()-1;
   Smaxj=seq2->length()-1;
   clearResult();
   //cerr << "global score: " << Smax << endl;
   return Smax;
}

template<class T>
int Dynaln<T>::local() {
   alntype=LOCAL;
   allocmem();
   int i,j,s1len, s2len;
   s2len=seq2->length();
   s1len=seq1->length();
   Smax=Smaxi=Smaxj=0;
//#ifdef DEBUG
   //cerr << "input sequences:\n" << *seq1 << *seq2 << endl;
   //cerr << "\ngap parameter inside local(): " << gapi << " " << gape << endl;
   //ST.show(cerr);
//#endif

   // fill the max array with zero
   // we only need the top row for memory
   // 0 as bourder value, 1 as first value, seq2len as last value
   for (j=0; j <= s2len; j++) M[j]=IX[j]=IY[j]=0;
   int leftM=0; 
   int leftIx=0;
   int leftIy=0;
   int CM=0, CIx=0, CIy=0, backptr, currScore; // current values
   /*
    * Use 3 arrays, plus left M,Ix,Iy to represent the 3 matrices: match, insert gap in X,
    * insert gap in Y.
    * X for column , Y row
    * o---- Y ---->
    * |
    * X
    * |
    *\|/
    * Here I have a little bit of code duplication for performance.
    * Use CM, CIx, CIy to decide the trace back pointer direction.
    * If there is a tie, then use maxima of diag, top, and left
    * to decide which direction to point to. Even with this
    * two layers of inforamtion, there is still cases where the
    * scores are identical for all three of them. The pointer is
    * saved as an OR operation, the it is left to the trace-back
    * algorithm to decide which direction to pick.
    */
   int maxDiag, maxTop, maxLeft;
   for (i=0; i < s1len; i++) {
#ifdef DEBUG
      cerr << "row i " << i << " Residue:" << (*seq1)[i] << endl;
#endif
      for (j=0; j<s2len; j++) {
         //cout << "col j " << j << " " << (*seq2)[j];
#ifdef DEBUG
         //cout << " " << j << (char)toupper((*seq2)[j]);
         //cerr << " " << j << " R:" << toupper((char)((*seq2)[j]));
         cerr << " " << j << " R:" << (*seq2)[j] << ' ';
#endif
         // j is the top row previus cell, C2[j] is the current sequence
         maxDiag = max(M[j], max(IX[j], IY[j]));
         CM=maxDiag + ST.lookup(C1[i], C2[j]);
         CIx=max(M[j+1]+gapi, IX[j+1]+gape); // above cell
         CIy=max(leftM+gapi, leftIy+gape);
         if (CM <= 0 && CIx <= 0 && CIy <= 0) {
            currScore = 0;
            //backptr = PTRNULL;
         }
         else if (CM > CIx && CM > CIy) { 
            currScore = CM;
            backptr = PTRDIAG;
         }
         else if (CIx > CM && CIx > CIy) {
            currScore = CIx;
            backptr = PTRTOP; // insert gap in seq2
         }
         else if (CIy > CIx && CIy > CM) {
            currScore = CIy;
            backptr = PTRLEFT; // insert gap in seq1
         }
         else if (CM == CIx && CM > CIy) {
            currScore = CIx;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            if (maxDiag > maxTop) backptr = PTRDIAG;
            else if (maxTop > maxDiag) backptr = PTRTOP;
            else {
               //cerr << i << "," << j 
               //   << " Top and Diag trace pointer have the same probability!\n";
               backptr = (PTRTOP | PTRDIAG);
            }
         }
         else if (CM == CIy && CM > CIx) {
            currScore = CIy;
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxDiag > maxLeft) backptr = PTRDIAG;
            else if (maxLeft > maxDiag) backptr = PTRLEFT;
            else {
               //cerr << i << ", " << j 
               //   << " Left and Diag trace pointer have the same probability!\n";
               backptr = (PTRLEFT|PTRDIAG); // this is very rare
            }
         }
         else if (CIy == CIx && CIy > CM) {
            currScore = CIy;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxTop > maxLeft) backptr = PTRTOP;
            else if (maxLeft > maxTop) backptr = PTRLEFT;
            else {
               //cerr << i << ", " << j
               //   << " Left and Top track pointer have the same probability!\n";
               backptr = (PTRLEFT | PTRTOP);
            }
         }
         else { // this state is rare but possible
            currScore = CM;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            /*
            cerr << i << ", " << j
               << " All three trace pointer have the same probability! very rare\n"
               << " Scores for CM, CIx, CIy: "
               << CM << ", " << CIx << ", " << CIy << endl
               << " Diag, Top, Left: " << maxDiag << ", " << maxTop << ", "
               << maxLeft << endl;
               */
            if (maxDiag > maxTop && maxDiag > maxLeft) backptr= PTRDIAG;
            else if (maxTop > maxDiag && maxTop > maxLeft) backptr= PTRTOP;
            else if (maxLeft > maxDiag && maxLeft > maxTop) backptr= PTRLEFT;
            else if (maxDiag == maxTop && maxDiag > maxLeft) {
               //cerr << "Diag and Top trace pointer have the same probability!\n";
               backptr= PTRDIAG | PTRTOP;
            }
            else if (maxDiag == maxLeft && maxDiag > maxTop) {
               //cerr << "Diag and Left trace pointer have the same probability!\n";
               backptr= PTRDIAG | PTRLEFT;
            }
            else if (maxTop == maxLeft && maxTop > maxDiag) {
               //cerr << "Top and Left trace pointer have the same probability!\n";
               backptr= PTRTOP | PTRLEFT;
            }
            else {
               //cerr << "it is very rare that the three values are identical: "
               //   << maxDiag << ", " << maxTop << ", " << maxLeft << endl;
               backptr= PTRDIAG | PTRTOP | PTRLEFT;
            }
         }
         if (currScore <= 0) backptr = PTRNULL;
#ifdef DEBUG
         cerr << " (" << CM << ", " << CIx << ", " << CIy << ") currScore=" << currScore;
#endif
         // save the trace-back pointer in the matrix S
         S[i*s2len+j] = backptr;
         if (currScore > Smax) {
            Smax=currScore;
            Smaxi=i; Smaxj=j;
         }
#ifdef DEBUG
         cerr << " max score: " << Smax << ", [" << Smaxi << ", " << Smaxj << "]" << endl;
#endif
         // ready for the next column 
         M[j]=leftM > 0? leftM : 0; 
         IX[j]=leftIx > 0? leftIx: 0; 
         IY[j]=leftIy > 0? leftIy : 0; // for the lower row
         leftM=CM > 0? CM : 0; 
         leftIx=CIx > 0? CIx : 0; 
         leftIy=CIy > 0? CIy : 0;
      }
      //cout << " max: " << Smax << " [" << Smaxi << ", " << Smaxj << "]\n";
      M[j]=CM > 0? CM : 0; 
      IX[j]=CIx > 0? CIx : 0; 
      IY[j]=CIy > 0? CIx : 0;
      leftM=leftIx=leftIy=0;
   }
   clearResult();
#ifdef DEBUG
   // show pointer matrix
   debug_showmatrix(cout);
   cerr << "max score " << Smax << endl;
#endif
   return Smax;
}


template<class T>
void Dynaln<T>::countnumgaps() {
   list<pair<int,int> >::const_iterator li=alnidx.begin();
   numgaps1=0;
   numgaps2=0;
   if (alnidx.empty()) return;
   while (li != alnidx.end()) {
      if (li->first == -1) {
         ++numgaps1;
         while (li != alnidx.end() && li->first == -1) ++li;
      }
      else ++li;
   }
   li=alnidx.begin();
   while (li != alnidx.end()) {
      if (li->second == -1) {
         ++numgaps2;
         while (li != alnidx.end() && li->second == -1) ++li;
      }
      else ++li;
   }
   //cerr << "gaps: " << numgaps1 << " " << numgaps2 << endl;
}

template<class T>
void Dynaln<T>::buildResult(const int delta1, const int delta2) {
   tracepointer(Smaxi, Smaxj);
   findAlnBoundary();
   buildAlnInfo(delta1, delta2);
   countnumgaps();
}

/*
 * marking at every 10 this could become a parameter if you
 * use very long sequence to do the alignment
 */
template<class T>
void Dynaln<T>::buildAlnInfo(const int delta1, const int delta2) {
   if (alnidx.empty()) {
      topaln.clear();
      bottomaln.clear();
      middle.clear();
      return;
   }
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i,j,counter;
   counter=0;
   char topchar, bottomchar;
   vector<int> topmark(alnidx.size(),-1);
   vector<int> bottommark(alnidx.size(),-1);

   while (lit != alnidx.end()) {
      i=lit->first; // top sequence index
      j=lit->second; // bottom sequence index
      //if (counter % 10 == 0) {
      if (counter % markevery == 0) {
         if (i>-1) { 
            topmark[counter]=i+1+delta1;
         }
         if (j>-1) {
            bottommark[counter]=j+1+delta2;
         }
      }
      if (i>=0) {
         topchar=toupper((*seq1)[i]); 
         //topchar=(*seq1)[i];
      }
      else {
         topchar=gapchar;
         ++gaplen1;
      }
      topaln += topchar;

      if (j>=0) {
         //bottomchar = (*seq2)[j];
         bottomchar = toupper((*seq2)[j]);
      }
      else {
         bottomchar = gapchar;
         ++gaplen2;
      }
      bottomaln += bottomchar;

      if (topchar != gapchar && bottomchar != gapchar) {
         if (topchar == bottomchar) {
            middle += idenchar;
            ++idencnt;
         }
         else if (ST.lookup(C1[i],C2[j]) >= simcut) {
            ++simcnt;
            middle += simchar;
         }
         else middle += ' ';
      }
      else {
         middle += ' ';
      }
      ++lit; ++counter;
   }
   string topr(alnidx.size(), ' ');
   string bottomr(alnidx.size(), ' ');
   string x;

   for (i=0; i<int(alnidx.size()); i++) {
      if (topmark[i] > -1) {
         x=itos(topmark[i]);
         topr.replace(i, x.length(), x);
      }
      if (bottommark[i] > -1) {
         x=itos(bottommark[i]);
         bottomr.replace(i, x.length(), x);
      }
   }
   topruler=topr;
   bottomruler=bottomr;
}

// the implementations are similar to getNgMatchArray
// I did not make a one depends on the other. The performance
// will be slower that way.
template<class T>
void Dynaln<T>::getMatchArray(vector<int> &marr, bool countsim) const {
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i,j,k,count;
   //vector<int> marr(alnidx.size(), 0);
   marr.resize(alnidx.size());
   k=0;
   while (lit != alnidx.end()) {
      i=lit->first;
      j=lit->second;
      if (i<0 || j<0) count= - 1;
      else if ((*seq1)[i] == (*seq2)[j]) count=2;
      else if (countsim && ST.lookup(C1[i],C2[j]) >= simcut) count=1;
      else count=0;
      marr[k++]=count;
      ++lit;
   }
}

// this version only cares about the non-gapped version
// this is being used for bootstrap
// this is used by bootstrap
template<class T>
void Dynaln<T>::getNgMatchArray(vector<int> &ngmarr, bool append, bool countsim ) const {
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i, j, count;
   if (!append) ngmarr.clear();
   while (lit != alnidx.end()) {
      i=lit->first;
      j=lit->second;
      if (i>=0 && j>=0) {
         if ((*seq1)[i] == (*seq2)[j]) count=2;
         else if (countsim && ST.lookup(C1[i],C2[j]) >= simcut) count=1;
         else count=0;
         ngmarr.push_back(count);
      }
      ++lit;
   }
}

template<class T>
ostream& Dynaln<T>::printSummaryHeader(ostream &ous, const string& dl, bool wantEntropy) {
   ous << "seq1name" << dl << "seq1length" << dl
      << "seq2name" << dl << "seq2length" << dl
      << "score" << dl << "identical" << dl << "similar" << dl
      << "alnlen" << dl << "numgaps1" << dl << "numgaps2" << dl
      << "gaplen1" << dl << "gaplen2" << dl
      << "seq1begin" << dl << "seq1end" << dl
      << "seq2begin" << dl << "seq2end";
   if (wantEntropy) {
      ous << dl << "seq1entropy" << dl << "seq2entropy";
   }
   return ous;
}

template<class T>
ostream& Dynaln<T>::printSummary(ostream &ous, const string &dl, bool wantEntropy) const {
   int ibase=1; // first base 1-based index
   ous << seq1->getName() << dl << seq1->length() << dl
      << seq2->getName() << dl << seq2->length() << dl
      << Smax << dl << idencnt << dl << simcnt << dl
      << alnidx.size() << dl
      << numgaps1 << dl << numgaps2 << dl
      << gaplen1 << dl << gaplen2 << dl
      << seq1begin+ibase << dl <<  seq1end+ibase << dl
      << seq2begin+ibase << dl <<  seq2end+ibase;
   if (wantEntropy) {
      ous << dl << setprecision(4) << (seq1->computeEntropy(seq1begin, seq1end)).first << dl
         << (seq2->computeEntropy(seq2begin, seq2end)).first << dl;
   }
   return ous;
}

template<class T>
void Dynaln<T>::printAlign(ostream &ous, const int w) const {
   //cerr << "print align info for human ...\n";
   ous << seq1->getName() << " x " << seq2->getName() << "  "
      << seq1begin << "-" << seq1end << "/" << seq1->length()
      << " | " 
      << seq2begin << "-" << seq2end << "/" << seq2->length()
      << endl
      << "Score=" << Smax 
      << " gap length: " << gaplen1 << " " << gaplen2
      << " num gaps: " << numgaps1 << " " << numgaps2
      << " idencnt=" << idencnt 
      << " simcnt=" << simcnt
      << " alnlen=" << topaln.length()
      << " identity=" << setprecision(5) << getIdentity()
      << " similarity=" << getSimilarity()
      << endl;
   // topaln, middle and bottomaln should be the same length
   if (topaln.length() != middle.length() 
         || middle.length() !=bottomaln.length()) {
      cerr << "the final alignment string is not right\n";
      exit(1);
   }
   //int i=0, i1=1, i2=1, charcnt;
   unsigned int i=0, j;
   //string tmp, tmp2, ruler1, ruler2;
   string ruler(w,' '), ruler1, ruler2;
   for (int x=0; x<w; x+=markevery) { // markevery default 10 residues
      ruler[x]='+';
   }
   while (i<topaln.length()) {
      j=i+w-1;
      while (isdigit(topruler[j])) ++j;
      ruler1=topruler.substr(i, j-i);
      j=i+w-1;
      while (isdigit(bottomruler[j])) ++j;
      ruler2=bottomruler.substr(i, j-i);

      ous << ruler1 << endl << ruler << endl
         << topaln.substr(i,w) << endl
         << middle.substr(i, w) << endl
         << bottomaln.substr(i,w) << endl << ruler << endl
         << ruler2 << endl << endl;
      i += w;
   }
}

template<class T>
string Dynaln<T>::alignToString(const int w) const {
   ostringstream ous;
   ous << seq1->getName() << " x " << seq2->getName() << "  "
      << seq1begin << "-" << seq1end << "/" << seq1->length()
      << " | " 
      << seq2begin << "-" << seq2end << "/" << seq2->length()
      << endl
      << "Score=" << Smax 
      << " gap length: " << gaplen1 << " " << gaplen2
      << " num gaps: " << numgaps1 << " " << numgaps2
      << " idencnt=" << idencnt 
      << " simcnt=" << simcnt
      << " alnlen=" << topaln.length()
      << " identity=" << setprecision(5) << getIdentity()
      << " similarity=" << getSimilarity()
      << endl;
   // topaln, middle and bottomaln should be the same length
   if (topaln.length() != middle.length() 
         || middle.length() !=bottomaln.length()) {
      cerr << "the final alignment string is not right\n";
      exit(1);
   }
   //int i=0, i1=1, i2=1, charcnt;
   unsigned int i=0, j;
   //string tmp, tmp2, ruler1, ruler2;
   string ruler(w,' '), ruler1, ruler2;
   for (int x=0; x<w; x+=markevery) { // markevery default 10 residues
      ruler[x]='+';
   }
   while (i<topaln.length()) {
      j=i+w-1;
      while (isdigit(topruler[j])) ++j;
      ruler1=topruler.substr(i, j-i);
      j=i+w-1;
      while (isdigit(bottomruler[j])) ++j;
      ruler2=bottomruler.substr(i, j-i);

      ous << ruler1 << endl << ruler << endl
         << topaln.substr(i,w) << endl
         << middle.substr(i, w) << endl
         << bottomaln.substr(i,w) << endl << ruler << endl
         << ruler2 << endl << endl;
      i += w;
   }
   return ous.str();
}


template<class T>
pair<string,int> Dynaln<T>::getNucleicConsensus() const {
   string tmp(bottomaln);
   int diff = 0;
   for (size_t i = 0; i<tmp.size(); ++i) {
      if (tmp[i] != topaln[i]) {
         if (tmp[i] == gapchar) tmp[i] = topaln[i];
         else {
            if ((tmp[i] == 'A' && topaln[i] == 'C') ||
                  (tmp[i] == 'C' && topaln[i] == 'A')) tmp[i] = 'M';
            else if ((tmp[i] == 'A' && topaln[i] == 'G')
                  || (tmp[i] == 'G' && topaln[i] == 'A')) tmp[i] = 'R';
            else if ((tmp[i] == 'A' && topaln[i] == 'T') ||
                  (tmp[i] == 'T' && topaln[i] == 'A')) tmp[i] = 'W';
            else if ((tmp[i] == 'C' && topaln[i] == 'G') ||
                  (tmp[i] == 'G' && topaln[i] == 'C')) tmp[i] = 'S';
            else if ((tmp[i] == 'C' && topaln[i] == 'T') ||
                  (tmp[i] == 'T' && topaln[i] == 'C')) tmp[i] = 'Y';
            else if ((tmp[i] == 'G' && topaln[i] == 'T') ||
                  (tmp[i] == 'T' && topaln[i] == 'G')) tmp[i] = 'K';
            ++diff;
         }
      }
   }
   return pair<string, int>(tmp, diff);
}

template<class T>
pair<string,string> Dynaln<T>::getNucleicConsensus2() const {
   string tmpt(topaln);
   string tmpb(bottomaln);
   int diff=0;
   for (size_t i = 0; i<tmpt.size(); ++i) {
      if (tmpt[i] == gapchar) tmpt[i] = tmpb[i];
      else if (tmpb[i] == gapchar) tmpb[i] = tmpt[i];
      else if (tmpt[i] != tmpb[i]) {
         ++diff;
      }
   }
   if (diff == 0) {
      return make_pair(tmpt, string());
   }
   else {
      return pair<string, string>(tmpt, tmpb);
   }
}

template<class T> vector<pair<char,int>> Dynaln<T>::getCigar1() const {
   vector<pair<char,int>> res;
   list<pair<int,int>>::const_iterator i=alnidx.cbegin();
   if (topBeginIndex() > 0) {
      res.push_back(make_pair('S', topBeginIndex()));
   }
   char C;
   int cnt;
   while (i != alnidx.cend()) {
      cnt=0;
      if (i->first > -1) {
         if (i->second > -1) {
            C='M';
            while (i != alnidx.cend() && i->first > -1 && i->second > -1) {
               ++cnt; ++i;
            }
         }
         else { // seq2 is gap
            C='D';
            while (i != alnidx.cend() && i->first > -1 && i->second == -1) {
               ++cnt; ++i;
            }
         }
      }
      else { // seq1 is gap, seq2 cannot be gap
         if (i->second == -1) {
            throw runtime_error("not possible for two gaps to algin");
         }
         C='I';
         while (i != alnidx.cend() && i->first == -1 && i->second > -1) {
            ++cnt; ++i;
         }
      }
      res.push_back(make_pair(C, cnt));
   }
   if (topEndIndex() < seq1->length()-1) { // end with soft clip
      res.push_back(make_pair('S', seq1->length()-1-topEndIndex()));
   }
   return res;
}

template<class T> 
vector<pair<char,int>> Dynaln<T>::getCigar2() const {
   vector<pair<char,int>> res;
   list<pair<int,int>>::const_iterator i=alnidx.cbegin();
   if (bottomBeginIndex() > 0) {
      res.push_back(make_pair('S', bottomBeginIndex()));
   }
   char C;
   int cnt;
   while (i !=  alnidx.cend()) {
      cnt=0;
      if (i->first > -1) {
         if (i->second > -1) {
            C='M';
            while (i != alnidx.cend() && i->first > -1 && i->second > -1) {
               ++cnt; ++i;
            }
         }
         else { // seq2 is gap
            C='I';
            while (i != alnidx.cend() && i->first > -1 && i->second == -1) {
               ++cnt; ++i;
            }
         }
      }
      else { // seq1 is gap, seq2 cannot be gap
         if (i->second == -1) {
            throw runtime_error("not possible for two gaps to algin");
         }
         C='D';
         while (i != alnidx.cend() && i->first == -1 && i->second > -1) {
            ++cnt; ++i;
         }
      }
      res.push_back(make_pair(C, cnt));
   }
   if (bottomEndIndex() < int(seq2->length()-1)) { // end with soft clip
      res.push_back(make_pair('S', seq2->length()-1-bottomEndIndex()));
   }
   return res;
}

template<class T>
void Dynaln<T>::clearResult() {
   // for global alignment if you have a special situation like 
   // ACCAAAAAACAACCCCCCCCCCCC------------------------
   // ------------------------GGGTTTTTTTGGGGGGGGGGGGTT
   // you will idencnt as zero, but you still need to clear the results.
   //if (idencnt == 0 && alntype == LOCAL) return;
   //The following version should apply to both local and global version.
   if (getAlnlen() == 0 ) return;
   idencnt=0;
   simcnt=0;
   numgaps1=0;
   numgaps2=0;
   gaplen1=0;
   gaplen2=0;
   //alnlen=0;
   topaln.clear();
   middle.clear();
   bottomaln.clear();
}

template<class T>
string Dynaln<T>::toDelimitedString(const string &dl, int ibase) const {
   ostringstream ous;
   try {
      /* production version
      ous << seq1->getName() << dl << seq1->length() << dl
         << seq2->getName() << dl << seq2->length() << dl
         << Smax << dl << idencnt << dl << simcnt << dl
         << alnidx.size() << dl
         << numgaps1 << dl << numgaps2 << dl
         << gaplen1 << dl << gaplen2 << dl
         << seq1begin << dl <<  seq1end << dl
         << seq2begin << dl <<  seq2end << dl
         << setprecision(4)
         << (seq1->computeEntropy(seq1begin, seq1end)).first << dl
         << (seq2->computeEntropy(seq2begin, seq2end)).first << dl;
      */
      // debug version, so we know which line crashed
      //pair<double,double> tmp=seq1->computeEntropy();
      //cerr << tmp.first << " " << tmp.second << endl;
      //cout << seq1begin << " " << seq1end << endl;
      ous << seq1->getName() << dl << seq1->length() << dl;
      ous   << seq2->getName() << dl << seq2->length() << dl;
      ous   << Smax << dl << idencnt << dl << simcnt << dl;
      ous   << alnidx.size() << dl;
      ous   << numgaps1 << dl << numgaps2 << dl;
      ous   << gaplen1 << dl << gaplen2 << dl;
      ous   << seq1begin+ibase << dl <<  seq1end+ibase << dl;
      ous   << seq2begin+ibase << dl <<  seq2end+ibase << dl;
      ous << setprecision(4) << (seq1->computeEntropy(seq1begin, seq1end)).first << dl;
      ous   << (seq2->computeEntropy(seq2begin, seq2end)).first << dl;
   }
   catch (exception &er) {
      cerr << er.what() << endl;
      cerr << "fetal error inside Dynaln::toDelimitedString()" << endl;
      exit(1);
   }
   return ous.str();
}

template<class T>
string Dynaln<T>::headers(const string &dl) {
   return "seq1name" + dl + "seq1len" + dl
      + "seq2name" + dl + "seq2len" + dl
      + "score" + dl
      + "identical" + dl + "similar" + dl + "alnlen" + dl
      + "seq1numgap" + dl + "seq2numgap" + dl
      + "seq1gaplen" + dl + "seq2gaplen" + dl
      + "seq1begin" + dl + "seq1end" + dl
      + "seq2begin" + dl + "seq2end" + dl
      + "seq1entropy" + dl + "seq2entropy";
}

template<class T>
int Dynaln<T>::getSeq1Length() const {
   if (seq1 == nullptr) throw runtime_error("NULL pointer for aligner's seq1");
   return seq1->length();
}

template<class T>
int Dynaln<T>::getSeq2Length() const {
   if (seq2 == nullptr) throw runtime_error("NULL pointer for aligner's seq2");
   return seq2->length();
}

template<class T>
pair<char, char> Dynaln<T>::getResiduesByTopPosition(int tpos, int &bpos) {
   list<pair<int, int> >::const_iterator itr = alnidx.begin();
   while (itr != alnidx.end()) {
      if ((tpos - 1) == itr->first) {
         bpos = itr->second + 1;
         if (itr->second == -1)
            return make_pair((*seq1)[itr->first], gapchar);
         else
            return make_pair((*seq1)[itr->first], (*seq2)[itr->second]);
      }
      ++itr;
   }
   bpos = -1;
   return pair<char, char>(gapchar, gapchar);
}

template<class T>
string Dynaln<T>::getBottomAlnByTopPosition(int tpos1, int tpos2) const {
   if (tpos1 < 0 || tpos2 > int(seq1->length())-1) {
      throw out_of_range("top positon cannot be out of seq1");
   }
   list<pair<int, int> >::const_iterator itr = alnidx.begin();
   int b=0, e=0;
   while (itr != alnidx.end()) {
      if ((tpos1 - 1) == itr->first) {
         break;
      }
      ++itr; ++b;
   }
   e=b;
   while (itr != alnidx.end()) {
      if ((tpos2 - 1) == itr->first) {
         break;
      }
      ++itr; ++e;
   }
   return getBottomAln().substr(b, e-b+1);
}

////////////////// Linear Space Algorithms ///////////////////////////////

// return the pointers
//void Dynaln::allocmemLS() {
template<class T>
void LSDynaln<T>::allocmem() {
   int Nr = this->seq1->length()+1;
   int Nc = this->seq2->length()+1;
   if (this->M == 0) {
      this->Ssize=2*Nc;
      this->numcol=Nc;
      //S = new int[Nc];
      this->M=new int[this->Ssize];
      this->IX = new int[this->Ssize];
      this->IY = new int[this->Ssize];
      //Itop = new int[Nc];  // top insert Y array
      //SR = new int[Nc];
      MR=new int[this->Ssize];
      IXR = new int[this->Ssize];
      IYR = new int[this->Ssize];
   }
   else {
      if (this->Ssize < 2*Nc) {
         //delete[] S; 
         delete[] this->M; delete[] this->IX; delete[] this->IY;
         //delete[] SR; 
         delete[] MR; delete[] IXR; delete[] IYR;

         this->Ssize=2*Nc;
         this->numcol=Nc;

         //S = new int[Nc];
         this->M=new int[this->Ssize];
         this->IX = new int[this->Ssize];
         this->IY = new int[this->Ssize];
         //Itop = new int[Nc];  // top insert Y array
         //SR = new int[Nc];
         MR=new int[this->Ssize];
         IXR = new int[this->Ssize];
         IYR = new int[this->Ssize];
      }
   }
   this->C1=this->seq1->getcode();
   this->C2=this->seq2->getcode();
}

// helper debug function
template<class T>
void showMatrix(pair<int*,int*> &M, int w) {
   for (int i=0; i<w; i++) {
      cerr << M.first[i] << " ";
   }
   cerr << endl;
   for (int i=0; i<w; i++) {
      cerr << M.second[i] << " ";
   }
   cerr << endl;
}
      
/* the correct version, using 4 arrays
 * M, IX, IY, S, and 4 for reverse
 * Separating the two algorithms, backward and forward
 * return the start index of top and bottom rows 
 * respectively.
 */
template<class T>
pair<int,int> LSDynaln<T>::computeScoreFW(int b1, int e1, int b2, int e2) 
{
   int Nc=e2-b2+2; // number of the clumns of the scoring matrix
   //int Nr=e1-b1+2;
   int i,j;
   this->M[0]=this->IX[0]=this->IY[0]=0;
   this->M[1]=this->IX[1]=this->IY[1]=this->gapi;
   for (j=2; j<Nc; j++) {
      this->IY[j]=this->IX[j]=this->M[j]=this->M[j-1]+this->gape;
   }
   this->M[Nc]=this->IX[Nc]=this->IY[Nc]=this->gapi;
   // matrix is 1 more than the sequence length
   int *MDia, *MLeft, *IXTop, *IYDia, *IXLeft, *IYLeft;
   int toprow=0;
   int bottomrow=Nc;
 /*
#ifdef DEBUG
   cout << "Matrix forward from computeScoreFW " << b1 << "-" << e1
      << " x " << b2 << "-" << e2 << "\n";
#endif
*/
   for (i=b1; i<=e1; i++) {
      MDia=this->M+toprow;  // MTop = MDia+1, no need to keep
      IYDia=this->IY+toprow; 
      IXTop=this->IX+toprow+1; //IXDia=IX+toprow; not needs
      MLeft=this->M+bottomrow;
      IXLeft=this->IX+bottomrow;
      IYLeft=this->IY+bottomrow;
      /*
#ifdef DEBUG
      cout << "row " << i-b1 << endl << *MLeft << "," << *IXLeft << "," << *IYLeft << " | ";
#endif
*/
      for (j=b2; j<=e2; j++) {
         *(MLeft+1)=max(*MDia, max(*(IXTop-1), *IYDia)) + this->ST.lookup(this->C1[i], this->C2[j]);
         *(IXLeft+1)=max(*(MDia+1) + this->gapi, *IXTop+this->gape);
         *(IYLeft+1)=max(*MLeft + this->gapi, *IYLeft+this->gape);
         /*
#ifdef DEBUG
         cout << ST.lookup(C1[i],C2[j]) << " " << (*seq1)[i] << "x" << (*seq2)[j] << " "
            << *(MLeft+1) << "," << *(IXLeft+1) << "," << *(IYLeft+1) << " | ";
#endif
*/
         ++MDia; ++IXTop; ++IYDia;
         ++MLeft; ++IXLeft; ++IYLeft;
      }
      /*
#ifdef DEBUG
      cout << endl;
#endif
*/
      // initialize the first column, boundary condition
      // This overwrites the top row first column if done
      this->M[toprow]=this->IX[toprow]=this->IY[toprow]=this->M[bottomrow]+this->gape;
      toprow=Nc-toprow;
      bottomrow=Nc-bottomrow;
   }
   // reverse the first column initialization effect
   // now the top and bottom has also been switched
   this->M[bottomrow]=this->IX[bottomrow]=this->IY[bottomrow]=this->M[toprow] - this->gape;
   //cout << "====================================\n";
   //return make_pair(toprow, bottomrow);
   //need to reverse the last operation
   return make_pair(bottomrow, toprow);
}

// (b1>e1 && b2>=e2) || (b1>=e1 && b2>e2)
template<class T>
pair<int,int> LSDynaln<T>::computeScoreBW(int b1, int e1, int b2, int e2) 
{
   int Nc=b2-e2+2;
   //int Nr=b1-e1+2;
   int i,j;
   MR[0]=IXR[0]=IYR[0]=0;
   MR[1]=IXR[1]=IYR[1]= this->gapi;
   for (j=2; j<Nc; j++) {
      IYR[j]=IXR[j]=MR[j]=MR[j-1]+this->gape;
   }
   MR[Nc]=IXR[Nc]=IYR[Nc]=this->gapi;
   // matrix is 1 more than the sequence length
   int *MDia, *MLeft, *IXTop, *IYDia, *IXLeft, *IYLeft;
   //int toprow, bottomrow, r;
   int toprow=0;
   int bottomrow=Nc;
   for (i=b1; i>=e1; i--) {
      MDia=MR+toprow;  // MTop = MDia+1, no need to keep
      IXTop=IXR+toprow+1; //IXDia=IX+toprow; not needs
      IYDia=IYR+toprow;
      MLeft=MR+bottomrow;
      IXLeft=IXR+bottomrow;
      IYLeft=IYR+bottomrow;
      //cout << "row " << b1-i << endl << *MLeft << ", " << *IXLeft << ", " << *IYLeft << " | ";
      for (j=b2; j>=e2; j--) {
         *(MLeft+1)=max(*MDia, max(*(IXTop-1), *IYDia)) + this->ST.lookup(this->C1[i], this->C2[j]);
         *(IXLeft+1)=max(*(MDia+1)+this->gapi, *IXTop+this->gape);
         *(IYLeft+1)=max(*MLeft+this->gapi, *IYLeft+this->gape);
         //copying and  ready for the next round
         ++MDia; ++IXTop; ++IYDia;
         ++MLeft; ++IXLeft; ++IYLeft;
      }
      //cout << endl;
      // initialize the first column, boundary condition
      MR[toprow]=IXR[toprow]=IYR[toprow]=MR[bottomrow]+this->gape;
      toprow=Nc-toprow;
      bottomrow=Nc-bottomrow;
   }
   // reverse the intialization effect (used for looping setup) when done
   MR[bottomrow]=IXR[bottomrow]=IYR[bottomrow]=MR[toprow]-this->gape;
   return make_pair(bottomrow, toprow);
}

template<class T>
int LSDynaln<T>::global() {
   if (this->seq1->length() < 1) {
      throw AlnInputException("seq1 empty");
   }
   else if (this->seq2->length() < 1) {
      throw AlnInputException("seq1 empty");
   }

   cerr << "\n*** Global linear space dynamic alignment ***\n";
   this->alntype=GLOBAL;
   allocmem();
   if (!this->alnidx.empty()) {
      this->alnidx.clear();
   }
   this->Smax=path(0,this->seq1->length()-1, 0, this->seq2->length()-1);
   this->clearResult();
}

// b1 < e1 && b2 < e2
template<class T>
int LSDynaln<T>::path(int b1, int e1, int b2, int e2) {
   int row=e1-b1+1; // length of seq1 range 
   int col=e2-b2+1;  // length of seq range
   pair<int, int> forward, backward;
   list<pair<int,int> > pathForward, pathBackward;
   list<pair<int,int> >::const_iterator it;
   int s1begin, s1end, s2begin, s2end; // divide boundary
   int scoreFB;
   //int i1,j1, i2,j2;
#ifdef DEBUG
   cout << "\nfinding best path for  " << b1 << "-" << e1
      << " x " << b2 << "-" << e2 << endl;
      //<< seq1->substr(b1+1, e1+1) << endl
      //<< seq2->substr(b2+1, e2+1) << endl;
#endif

   if (row < 2 && col < 2) {
      //cout << b1 << "," << b2 << " |\n";
      this->alnidx.push_back(make_pair(b1,b2));
      return this->ST.lookup((*(this->seq1))[b1], (*(this->seq2))[b2]);
   }
   else if (row < 2 || col < 2) {
      forward=computeScoreFW(b1,e1,b2,e2);
      s1end=e1;
      s2end=e2;
      pathForward=tracebackFW(forward,b1,s1end,b2,s2end);
      this->alnidx.insert(this->alnidx.end(), pathForward.begin(), pathForward.end());
      int bb=forward.second;
      return max(this->M[bb+col], max(this->IX[bb+col], this->IY[bb+col]));
   }

   s1end=(b1+e1)/2;  // top part end
   s1begin=s1end+1; // bottom part begin
   forward=computeScoreFW(b1, s1end, b2, e2);
   backward=computeScoreBW(e1, s1begin, e2, b2);
   //int *F=S+forward.second;
   //int *B=S+backward.second;
   ///////// forward ///////////
   int *scoreM=this->M+forward.second;
   int *scoreIX=this->IX+forward.second;
   int *scoreIY=this->IY+forward.second;
   ////// backward //////////////
   int *scoreMR=MR+backward.second;
   int *scoreIXR=IXR+backward.second;
   int *scoreIYR=IYR+backward.second;

   int maxScore=-99999999; // a very small number
   int maxj, j;
   maxj=-1;
   //cout << "find the division cell\n";
   for (j=0; j<=col; j++) {
      scoreFB=max(scoreM[j], max(scoreIX[j], scoreIY[j]))
            + max(scoreMR[col-j], max(scoreIXR[col-j], scoreIYR[col-j]));
      /*
#ifdef DEBUG
      cout << scoreFB << " " << max(scoreM[j], max(scoreIX[j], scoreIY[j]))
            << " + " << max(scoreMR[col-j], max(scoreIXR[col-j], scoreIYR[col-j]))
            << " | ";
#endif
*/
      //if (F[j] + B[col-j] > maxScore) {
      if (scoreFB > maxScore) {
         maxScore= scoreFB;
         maxj=j;
      }
   }
   //cout << endl;
   // divide the large problem into the following two smaller ones
   // The index is in terms of the sequences
   // seq1 b1----s1end  s1begin----e1
   //      |||||||||||  |||||||||||||
   // seq2 b2----s2end  s2begin----e2
   // The traceback functions works on indices on the matrix which
   // has an additional dimention
   /*
   if (maxj == 0) {
      cerr << "Abnormal division result with maxscore: " << maxScore
         << " and maxj: " << maxj << endl
         << seq1->substr(b1+1, e1+1) << "\n" << seq2->substr(b2+1,e2+1)
         << endl << "seq1 divisions " << s1end << " | " << s1begin << endl;
      debug_showmatrixForward(forward,col+1);
      debug_showmatrixBackward(backward,col+1);
      exit(1);
   }
   */
   s2end=b2+maxj-1;
   s2begin=e2-col+maxj+1;
#ifdef DEBUG
   cout << "new seq1 division point: " << s1end << " " << s1begin << endl;
   cout << "new seq2 division point: " << s2end << " " << s2begin << endl;
   // speed up by skipping identical residues
#endif
   bool topDone=false;
   if (maxj == 0)  {
      topDone=true;
      while (s1end >= b1) {
         pathForward.push_front(make_pair(s1end--, -1));
      }
   }
   else if (s2end == b2 && b1 == s1end) {
      topDone=true;
      pathForward.push_front(make_pair(b1,b2));
   }
   else {
      //cout << "Before tracebackFW s1end,s2end(" << s1end << ", " << s2end << ")\n";
      pathForward = tracebackFW(forward, b1, s1end, b2, s2end);
      //cout << "After tracebackFW s1end,s2end(" << s1end << ", " << s2end << ")\n";
      if (s1end <= b1 || s2end <= b2) {
         //cout << "the new ends indicates finish of job after tracebackFW! "
            //<< b1 << "--" << s1end << " |x| " 
            //<< b2 << "--" << s2end << "\n";
         topDone=true;
         //exit(1);
      }
   }
#ifdef DEBUG
   cerr << "the ends before tracebackBW! "
      << s1begin << "--" << e1 << " || " 
      << s2begin << "--" << e2 << "\n";
   if (s2begin > e2) {
      cerr << "Division finished the bottom problem with maxscore: " << maxScore
         << " and maxj: " << maxj << endl
         << this->seq1->substr(b1+1, e1+1) << "\n" << this->seq2->substr(b2+1,e2+1)
         << endl << "seq1 divisions " << s1end << " | " << s1begin << endl;
      debug_showmatrixForward(forward,col+1);
      debug_showmatrixBackward(backward,col+1);
   }
#endif
   bool bottomDone=false;
   if (maxj == col) { // s2begin == e2+1
      bottomDone=true;
      while (s1begin <= e1) {
         pathBackward.push_back(make_pair(s1begin++, -1));
      }
   }
   else if (s2begin == e2 && e1==s1begin) {
      bottomDone=true;
      pathBackward.push_back(make_pair(e1,e2));
   }
   else {
      pathBackward = tracebackBW(backward, e1, s1begin, e2, s2begin);
      if (s1begin >= e1 || s2begin >= e2) {
         //cerr << "new ends indicated finish of the bottom job after tracebackBW! "
            //<< s1begin << "--" << e1 << " || " 
            //<< s2begin << "--" << e2 << "\n";
         bottomDone=true;
         //exit(1);
      }
   }
   pathForward.insert(pathForward.end(), pathBackward.begin(), pathBackward.end());
   if (!topDone) path(b1, s1end, b2, s2end);
   // save the result into the shared alnidx list
   this->alnidx.insert(this->alnidx.end(), pathForward.begin(), pathForward.end());
   if (!bottomDone) path(s1begin, e1, s2begin, e2);
   return maxScore;
}

/* forward version, b1<e1 && b2<e2
 * the matrix has one additional row,column than the sequences
 */
template<class T>
list<pair<int,int> > LSDynaln<T>::tracebackFW(pair<int, int> mrow, int b1, int &e1, int b2, int &e2) 
{
   list<pair<int, int> > path;
   int *MTop=this->M+mrow.first;
   int *MButtom=this->M+mrow.second; 
   int *IXTop=this->IX+mrow.first; 
   int *IXButtom=this->IX+mrow.second;
   int *IYTop=this->IY+mrow.first;
   int *IYButtom=this->IY+mrow.second;
   int r,c;
   r=e1-b1+1;
   c=e2-b2+1; // column index in the matrix

   int R=r;  // remember the last row
   //int i,j; // index in sequence for reverse lookup
   // stop the loop after moving into the top row
   //i=e1;
   //j=e2;
   int score;
   // stop when moved to the top row.
   while (c>0 && R-r != 1 && e2 >= b2) {
      score=max(MButtom[c], max(IXButtom[c], IYButtom[c]));
      if (score == MButtom[c]) {
         path.push_front(make_pair(e1, e2));
         --r; --c;
         --e1; --e2;
      }
      else if (score == IXButtom[c]) {
         path.push_front(make_pair(e1, -1));
         --r; --e1;
      }
      else if (score == IYButtom[c]) { 
         // from left, trace till it ends
         path.push_front(make_pair(-1, e2));
         --c; --e2;
      }
      else {
         cerr << "track back error. not possible\n";
         exit(1);
      }
   }
   // we are not working in the top row, and only move to the left
   while (c>0 && R-r == 1 && e2>=b2) {
      score=max(MTop[c], max(IXTop[c], IYTop[c]));
      if (score==MTop[c]) {
         path.push_front(make_pair(e1, e2));
         --r; --c; --e1; --e2;
      }
      else if (score == IXTop[c]) {
         path.push_front(make_pair(e1, -1));
         --r; --e1;
      }
      else if (score == IYTop[c]) {
         path.push_front(make_pair(-1, e2));
         --c; --e2;
      }
   }
   while (c==0 && r>0 && e1 >= b1) {
      path.push_front(make_pair(e1--, -1));
      --r;
   }
   //if (i>=b1) e1=i;
   //if (j>=b2) e2=j;
   return path;
}

// backward version, b1 >= e1 && b2 >= e2
// at this point having two separate version is better than keeping one for both
// Essentially copying and pasting
template<class T>
list<pair<int,int> > LSDynaln<T>::tracebackBW(pair<int, int> mrow, int b1, int &e1, int b2, int &e2) 
{
   list<pair<int, int> > path;
   int *MTop=this->M+mrow.first;
   int *MButtom=this->M+mrow.second; 
   int *IXTop=this->IX+mrow.first; 
   int *IXButtom=this->IX+mrow.second;
   int *IYTop=this->IY+mrow.first;
   int *IYButtom=this->IY+mrow.second;
   int r,c;
   r=b1-e1+1;
   c=b2-e2+1;

   int R=r;  // remember the last row
   //int i,j;  // index in sequence for reverse lookup
   // stop the loop after moving into the top row
   //i=e1;
   //j=e2;
   int score;
   // stop when moved to the top row.
   while (c>0 && R-r < 1 && e2 <= b2) {
      score=max(MButtom[c], max(IXButtom[c], IYButtom[c]));
      if (score == MButtom[c]) {
         path.push_back(make_pair(e1, e2));
         --r; --c; ++e1; ++e2;
      }
      else if (score == IXButtom[c]) {
         path.push_back(make_pair(e1, -1));
         --r; ++e1;
      }
      else if (score == IYButtom[c]) { 
         // from left, trace till it ends
         path.push_back(make_pair(-1, e2));
         --c; ++e2;
      }
      else {
         cerr << "track back error. not possible\n";
         exit(1);
      }
   }
   // we are not working in the top row, and only move to the left
   while (c>0 && R-r == 1 && e2 <= b2) {
      score=max(MTop[c], max(IXTop[c], IYTop[c]));
      if (score==MTop[c]) {
         path.push_back(make_pair(e1, e2));
         --r; --c; ++e1; ++e2;
      }
      else if (score == IXTop[c]) {
         path.push_back(make_pair(e1, -1));
         --r; ++e1;
      }
      else if (score == IYTop[c]) {
         path.push_back(make_pair(-1, e2));
         --c; ++e2;
      }
   }
   // special case trace when column 0 of the matrix is reached
   while (c==0 && r>0 && e1 <= b1) {
      path.push_back(make_pair(e1++, -1));
      --r;
   }
   //if (i<=b1) e1=i;
   //if (j<=b2) e2=j;
   return path;
}

// overwrite the square version
template<class T>
void LSDynaln<T>::buildResult(const int delta1, const int delta2) {
   //cerr << " *** calling LSDynaln buildResult() \n";
   list<pair<int, int> >::const_iterator lit= this->alnidx.begin();
   int i,j;
   char topchar, bottomchar;
   //clearResult();
   // find sequence begin and end, removing intial gaps
   // The local version produce this result in the
   // matrix computation step.
   if (this->alntype==GLOBAL) {
      this->seq1begin=this->seq2begin=0;
      this->seq1end=this->seq1->length()-1;
      this->seq2end=this->seq2->length()-1;
      if (lit->first == -1) {
         while (lit != this->alnidx.end() && lit->first == -1) ++lit;
         this->seq1begin=lit->first; // should be 0
         this->seq2begin=lit->second;
         if (this->seq2begin == -1) {
            cerr << "Inside buildResult(), indel of the two aligned sequence follow each others, the algorithme might not be perfet for the linear space algorithm!\n";
            list<pair<int,int> >::const_iterator lix=lit;
            while (lix != this->alnidx.end() && lix->second == -1) {
               ++lix;
            }
            this->seq2begin=lix->second;
         }
      }
      else if (lit->second == -1) {
         while (lit != this->alnidx.end() && lit->second == -1) ++lit;
         this->seq1begin=lit->first;
         this->seq2begin=lit->second;  // should be zero
         if (this->seq1begin == -1) {
            cerr << "Inside buildResult(), indel of the two aligned sequence follow each others, the algorithme might not be perfet for the linear space algorithm!\nSeq1 start with -1\n";

            list<pair<int,int> >::const_iterator lix=lit;
            while (lix != this->alnidx.end() && lix->first == -1) {
               cerr << lix->first << "x" << lix->second << "  ";
               ++lix;
            }
            this->seq1begin=lix->first;
            cerr << "\nthe new seq1begin: " << this->seq1begin << endl;
         }
      }

      list<pair<int,int> >::reverse_iterator ri=this->alnidx.rbegin();
      if (ri->first == -1) {
         while (ri != this->alnidx.rend() && ri->first == -1) ++ri;
         this->seq1end=ri->first;
         this->seq2end=ri->second;
      }
      else if (ri->second == -1) {
         while (ri != this->alnidx.rend() && ri->second == -1) ++ri;
         this->seq1end=ri->first;
         this->seq2end=ri->second;
      }
   }
   this->buildAlnInfo(delta1, delta2);
   this->countnumgaps();
}

/* use two rows of array to memorize the coordinate
 * of the alignment starting point
 */
template<class T>
int LSDynaln<T>::local() {
   this->alntype=LOCAL;
   allocmem();
   const int col=this->seq2->length()+1; // column of matrix
   const int row=this->seq1->length()+1;
   int i,j;
   // initialize match and indel score for row 0
   for (j=0; j< col; j++) {
      this->M[j]=this->IX[j]=this->IY[j]=0;
      // initialize position to default value -1
      MR[j<<1]=-1;  // faster version
      MR[(j<<1)+1]=-1;
   }
   this->M[col]=this->IX[col]=this->IY[col]=0;
   int *Mdia, *Mleft, *Mtop;
   int *IXdia, *IXleft, *IXtop;
   int *IYdia, *IYleft, *IYtop;
   // use MR array to remember the position
   // in the scoring matrix, the first column was set to -1

   int PLi, PLj, PCi, PCj; // position left and current
   //int *It;  // pointer to top insertion
   this->Smax=0;
   // for tracing the starting point of the alignment
   int *pos; // pointer to the position
   int Pmaxi, Pmaxj, score, scoreM, scoreIX, scoreIY;
   int top=0;
   int bottom=col;

   for (i=0; i < row-1; i++) {
      Mdia=this->M+top; Mleft=this->M+bottom;
      IXdia=this->IX+top; IXleft=this->IX+bottom;
      IYdia=this->IY+top; IYleft=this->IY+bottom;
      PLi=-1; PLj=-1;
      pos=MR+2;
      for (j=0; j < col-1; j++) {
         *(Mleft+1)=scoreM=max(*Mdia, max(*IXdia, *IYdia)) + this->ST.lookup(this->C1[i], this->C2[j]);
         *(IXleft+1)=scoreIX=max(*(Mdia+1)+this->gapi, *(IXdia+1)+this->gape);
         *(IYleft+1)=scoreIY=max(*Mleft+this->gapi, *IYleft+this->gape);
         score=max(max(scoreM, scoreIX), max(0,scoreIY));
         if (score > 0) {
            // not looking back on cell, the start position
            // is one before the actual alignment.
            if (score == scoreM) {
               PCi=*(pos-2); PCj=*(pos-1); // diagnal
            }
            else if (score == scoreIX) {
               PCi=*pos; PCj=*(pos+1);  // top
            }
            else if (score == scoreIY) {
               PCi=PLi; PCj=PLj; // left
            }
            else {
               cerr << "wrong condition in memorizing the starting point of local alignment\n";
               exit(1);
            }
            if (score > this->Smax) {
               this->Smax = score;
               this->Smaxi=i; this->Smaxj=j;
               Pmaxi=PCi; Pmaxj=PCj;
            }
         }
         else {
            PCi=i; PCj=j;
         }
         *(pos-2)=PLi; *(pos-1)=PLj; // left overwrite diag 
         PLi=PCi; PLj=PCj;       // current overwrite left
         pos += 2; // move two spaces to the right
         ++Mdia; ++Mleft;
         ++IXdia; ++IXleft;
         ++IYdia; ++IYleft;
      }
      this->M[top]=this->IX[top]=this->IY[top]=this->M[bottom]+this->gape;
      top=col-top;    // mechanism to reverse between 0 and 1
      bottom=col-bottom;
   }
   this->seq1begin=Pmaxi+1;
   this->seq1end=this->Smaxi;
   this->seq2begin=Pmaxj+1;
   this->seq2end=this->Smaxj;
   if (!this->alnidx.empty()) this->alnidx.clear();
   path(this->seq1begin, this->seq1end, this->seq2begin, this->seq2end);
   this->clearResult();
   return this->Smax;
}

template<class T>
LSDynaln<T>::~LSDynaln<T>() {
   //delete[] S; 
   //delete[] Itop;
   //delete[] M;
   //delete[] IX;
   //delete[] IY;
   // The base class members are done by base class
   delete[] MR;
   delete[] IXR;
   delete[] IYR;
}

// take top and bottom row pointers
template<class T>
void LSDynaln<T>::debug_showmatrixForward(pair<int,int> &rows, const int col) {
   int *Matchptr, *IXptr, *IYptr;
   Matchptr = this->M + rows.first;
   IXptr = this->IX + rows.first;
   IYptr = this->IY + rows.first;
   int j;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
   Matchptr = this->M + rows.second;
   IXptr = this->IX + rows.second;
   IYptr = this->IY + rows.second;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
}
template<class T>
void LSDynaln<T>::debug_showmatrixBackward(pair<int,int> &rows, const int col) {
   int *Matchptr, *IXptr, *IYptr;
   Matchptr = MR + rows.first;
   IXptr = IXR + rows.first;
   IYptr = IYR + rows.first;
   int j;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
   Matchptr = MR + rows.second;
   IXptr = IXR + rows.second;
   IYptr = IYR + rows.second;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
}
}
#endif
