#ifndef BIOSEQ_H
#define BIOSEQ_H

// (C) 2002 Kemin Zhou at orpara.com

/** sequence is stored as string if needed can be converted to
 * different formats
 * Sequence id is also converted to string, even it could be
 * integer id.  This is because humans usually like to use
 * strings. But in the future this could be changed to integers.
 * It will speed up computation.
 */
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <list>
#include <vector>
#include <range.h>
#include <interval.h>

#include "codon.h"

using namespace std;

namespace orpara {
class greaterRangeLength {
   public:
      bool operator()(const Range &r1, const Range &r2) const {
         return r1.length() > r2.length();
      }
};

/** This exception type is define to be used
 * by the bioseq objects
 */
class bioseqexception : public exception {
   public:
      bioseqexception() throw() : errstr("bioseq error") { }
      bioseqexception(const string &msg) throw() : errstr(msg) { }
      ~bioseqexception() throw() { }
      const char* what() const throw() { return errstr.c_str(); }
   private:
      string errstr;
};

// some helper functions
/** class wide function that can be used without constructing
 * a bioseq object.  This will increase the performance if 
 * you are only interested in the operation but not in using
 * the bioseq object and its derived classes
 */
void printFasta(ostream &ous, const string& seq, unsigned int width=70);
/** operation only make sense on DNA or RNA sequences 
 * reverse complement input string seq
 * @param seq the input sequence.
 * @return the reverse complement of seq
 */
string reverseComplement(const string& seq);
/** 
 * a version for better performance 
 * Change the sequence.
 * @param seq the input. After this operation seq will
 *    be modified to become the reverse-complement.
 */
void reverseComplementInPlace(string& seq);
/** 
 * load all the sequences in the file into 
 * a map for later usage. Only the id is used,
 * title information is not cached.
 *
 * @return true for success, false for failure.
 */
bool loadFastaIntoMap(const string &file, map<string,string> &store);
/** This file also defines a sequence type
 * as enum so we can save storage space
 */
enum SequenceType {PROTEINSEQ=1, DNASEQ=2, RNASEQ=4, NUCLEIC_ACID=6, GENERIC=7, DNASEQQUAL=8, DNASEQQUALCOUNT=9, UNKNOWN=0};
/** Algorithms about ORF:
  * prefix ORF ORF starting from 5'-END without start codon.
  * suffix ORF ORF end at 3' end without stop codon
  * full ORF ORF with start and stop codons
  * Protein index to RNA index transformation:
  * in 0-based index system
  *  Ri=frame + 3*Pi First base of codon
  * For end of stop codon, you need to add 2
  */
/** helper function to be used by other methods
 * Uses more basic type, for convinience
 * use 1-based index inclusive.
 * @param begin start position of translation, first base.
 * @param end end position of translation, third base of codon if
 *     not at the end of the sequence.  This can be any of the
 *     codon bases. It just generate partial peptides with the
 *     last amino acid unspecified.
 *     end of zero signify the end of the sequence (All the way to the end).
 */
void translate(string &pep, const string &seq, unsigned int begin, unsigned int end=0);
//void translate(string &pep, const string &seq, int begin, int end=0);
//   translate(pep, seq, (unsigned)begin, (unsigned)end); }
/** global function to be used in a 
 * more flexible way
 */
int countInternalStops(const string &seq);
/** find all ORF 
 * (in the middle of pepseq M...*) or ...* or M...
 * in all three reading frames of rna, and pick the 
 * longest one. Set the pep seq, and b, e as 0-based index
 * inclusive [b,e] in RNA index
 * @param f is the frame [0,1,2]
 * actually frame can be derived from b % 3
 * frame is always b%3, actual begin is 0 if b < 3 
 * and start is not M
 */
void longestORFPlus(const string &rna, string &pep, int &b, int &e, int &f);
void longestORFPlus(const string &rna, string &pep, int &b, int &e);
/** Suffix ORF is ORF with start but not stop */
void longestORFPlusSuffix(const string &rna, string &pep, int &b, int &e);
/** Prefix ORF is the one with stop but no start */
void longestORFPlusPrefix(const string &rna, string &pep, int &b, int &e);
/** Find the longest ORF of
 * 1----*, and M---*, return both of them in one operation
 * Let the user decide what to do with the two values.
 * The * stands for stop codon. 
 * use 0-based index, inclusive [b,e]
 * e is the third Base of the stop codon.
 * b is the first base of the start codon if complete ORF. If nostart
 * orf, the b is the frame, implying start from 0.
 * b,e is packed into the pair data structure.
 *
 * nterm, and full contain the [b,e] in RNA coordinates.
 *
 * If the sequence does not contain NoStart or Full ORF, then return
 * false;
 * @return true if contain at least one of NoStart or Full ORF
 *         false if both are missing. In this situation, the 
 *         pep will be set to the no-start-and-no-stop ORF is
 *         it exists and full will be set to NostartNostrop frame
 *         to the end of the rna seq. This usually happens when 
 *         the sequence is short.
 */
bool longestNoStartORFPlus(const string &rna, pair<int,int> &nterm, string &npep,
      pair<int,int> &full, string &pep);
bool longestNoStopORFPlus(const string &rna, pair<int,int> &cterm, string &cpep,
      pair<int,int> &full, string &pep);


/** This translates full ORF and prefix ORF, for suffix ORF
  * you need special treatment, from protein space
  * to RNA space. For prefix, will use the start to encode frame
  * information; if start < 3 then it contain frame
  * info. Full ORF starts with M.
  * @param frame [0,1,2]
  */
Range P2Rindex(const Range &ofb, const int frame);
Interval P2Rindex(const Interval &ofb, const int frame);
int P2Rindex(const int pos, const int frame);
//{ return 3*pos+frame; }

/** helper used by findAllORFIndex */
bool overlayRVR(const Range &r, const vector<Range> &vr);
/** @return non-overlapping ORFs that are supposed to be more 
  * biological. some of the ORF might not be the longest,
  * but the overall Coding is optimal. The algorithm is simply
  * a packing algorithm by selecting the largest first,
  * then select the next non-overlapping ORF.
  * You need to do your own filtering.
  * The returned range are sorted from small to large
  *  regardless of direction using the less operator
  *  of Range.
  * @param base [0,1] default 1
  *
  * The Range returned is 0-based index.
  * this method should be applied to sequence at least 
  * 100 nt long.
  */
vector<Range> findAllORFIndex(const string &rna, int base=1,
      int pep_lencut=22, int HHcut=160, int HTcut=120, int TTcut=10);
/** find all the ORF index bounds in 0-based index
  * Number is in protein space
  * minimum AA length to register.
  * 25 aa is the minimum we are goint to register.
  * @param orfrange is the result of this operation. It will make
  *     it empty at the begining of the run if it was not.
  *     The output is in ss index. the end of the range is the 
  *     '*' symbol for complete ORF or prefix ORF. The begin
  *     of the range is the index of 'M'. 
  * @param ss is the input peptide sequence
  * @cuotff is a integer number. Peptide shorter than this are
  *     ignored in the searching phase. This parameter should
  *     not affect the real performance of this algorithm.
  *     It provided fine control. Currently I am using 25 aa.
  *     partial ORF on either end of the DNA are given
  *     1/2*cutoff or 35 aa limits whichever is larger.
  *     This is based on the average random ORF length of 21 aa.
 */
void findAllPepORFIndex(list<Range> &orfrange, const string &ss, unsigned int cutoff);

/** 
 * used by encode function for efficient lookup. 
 * hashbase function is declared in codon.h
 *
 * The aachar2num assumes a is in upper letter.
 * The algorithm is a - 'A'.
 * The bioseq used the same algorithm.
 * Cannot be inline, must be in implementation
 * file, otherwise get multiple definition error.
 */
int aachar2num(char a);

/** convert amino acid number to character */
char aanum2char(int c);


/** the bioseq base class describing the biological sequence.
 */
class bioseq {
   public:
      bioseq() : seq(), name(), title(), code(nullptr) {}
      /**
       * copy constructor.
       * Will not copy code.
       */
      bioseq(const bioseq &s) 
         : seq(s.seq), name(s.name),
            title(s.title), code(nullptr) { }
      /**
       * Move constructor from another sequence.
       */
      bioseq(bioseq &&s) 
         : seq(std::move(s.seq)), 
           name(std::move(s.name)), title(std::move(s.title)),
           code(s.code)
      { 
         if (s.code != nullptr) s.code=nullptr; 
         //if (!s.title.empty()) title=std::move(s.title); 
      }

      /**
       * Constructor from string
       */
      bioseq(const string &s) 
         : seq(s), name(), title(), code(nullptr) {}
      bioseq(string &&s) 
         : seq(std::move(s)), name(), title(), code(nullptr) {}
      /** 
       *  Constructor from name and string sequence
       *  @param n name of the bioseq
       *  @param s sequence 
       */
      bioseq(const string &n, const string &s) 
         : seq(s), name(n), title(), code(nullptr) {}
      bioseq(string &&n, string &&s) 
         : seq(std::move(s)), name(std::move(n)), title(), code(nullptr) {}
      /** 
       * @param n name
       * @param s sequence
       * @param t title
       */
      bioseq(const string &n, const string &s, const string &t) 
         : seq(s), name(n), title(t), code(nullptr) {}
      //virtual ~bioseq() { 
      virtual ~bioseq() { if (code != nullptr) delete[] code; }

      /**
       * restore this object to pure empty state
       */
      void clear() { 
         if (code != nullptr) {
            delete[] code; 
            code=nullptr; 
         }
         seq.clear(); name.clear(); title.clear(); 
      }

      /** 
       * change the underlying sequence
       */
      void assign(const string &s) { 
         seq=s; 
         if (code != nullptr) {
            delete[] code;
            code=nullptr;
         } 
      }
      void assign(string&& s);
      /**
       * alias for assign
       * Update the underlying sequence with s.
       * @param s input sequence
       */
      void setSequence(const string &s) { assign(s); }
      void setSequence(string &&s);

      ///////////////// input function //////////////////////
      
      /** 
       * Read only one sequence from fasta formated files.
       * Should not do repeated calling of this function.
       * It will read only the first sequence even if you have 
       * multiple sequences in the file.
       **/
      bool read(const string &file);

      /** 
       * Repeated call of this method is fine, and this object will become the
       * next sequence stored in the underlying file.
       * @param header you need to provide a header to let the biosequence to
       *        rember the last header read. An empty string should be provided
       *        initially.
       * @return false if no more sequence in the input stream.
       *         true if still got more
       *
       * Note: for better design, there should be a BioseqReader class.
       *       I may add one later.
       */
      bool read(istream &ins, string &header);
      /**
       * Get one sequence with this call. Can do repeat call of 
       * this function in a loop.
       * @return true if not reaching the end of file.
       *    fasle if reached the end of file. Will throw
       *    exception if bad fasta file encountered.
       */
      bool read(istream &ins);
      /** 
       * Essentially the same as read. Just a different foramt.
       */
      friend istream& operator>>(istream& ins, bioseq &sq);

      /**
       * return a reference to the underlying sequence string
       * the input string fas is in fasta format.
       * removes all non residue character.
       */
      string& importFasta(const string &fas);

      //////// output functions ////////////////////////////

      /** 
       * makes a copy of the underlying sequence
       * and return the sequence as string
       *
       * @return a copy of the underlying sequence as string 
       *    that can be altered without affecting this object.
       * @see getSequence()
       */
      string toString() const { return seq; }

      /** 
       * use the name as identifier, no title information
       */
      void printFastaWithHeader(ostream &ous, unsigned int width=70) const;
      
      /** 
       * only print the sequence, no id or title.
       * Use printFastaWithHeader method for file output.
       * It will print the last line terminator.
       * */
      void printFasta(ostream &ous, unsigned int width=70) const;
      //void printFasta(ostream &ous, int width=70) {
      //   printFasta(ous, (unsigned)width); 
      //}
      /** 
       * fasta format, 70 residues per line.
       * If the sequence has name it will produce a full
       * fasta file, otherwise only the sequence.
       * Not virtual function, different derived class may not
       * be able to print more information.
       * The end of the fasta sequence is Terminated.
       * Caller don't need to output ENDL.
       */
      friend ostream& operator<<(ostream &ous, const bioseq &s);
      /** this is used for polymorphic behavior of the 
       * derived classes.
       */
      virtual ostream& print(ostream &ous, int width=80) const;

      /** The biologist way of thinking.
       * inclusive range [b,e], different from C++ substr!
       * 1-based index
       * @return a copy of the underlying substring.
       */
      string substr(int b, int e) const; 
      /**
       * @return the substring of the underlying string object
       *    starting from 1-based index b to the end.
       */
      string substr(int b) const; 
      /**
       * C-style index, same as C::substr() function
       * @param b 0-based index of start
       * @param len length of the string.
       * @return a copy of the underlying substring.
       */
      string substring(const int b, const int len) const; 
      /**
       * Same function as C++ string.substr
       * @return the substring from b to end of the underlying string object.
       * @param b start of the substring 0-based index
       */
      string substring(const int b) const; 

      /** performance is very bad if sequence is long.
       * 1-based index, inclusive [b, e]
       * If b > e this function will reverse complenent the sequence.
       * The name of the returned sequence will be derived from the original
       * sequence with range in it.
       */
      bioseq subseq(int b, int e=-1) const;

      /** use the C-style indexing
       * This is an expensive operation. It is making a copy
       * of the subsequence.
       *
       * For efficiency, the name of the sequence was not returned.
       * The subsequence will have empty name.
       * @param b 0-based begin index
       * @param len if not given assuming to the end of sequence.
       */
      bioseq subsequence(int b, int len=-1) const;
      /** use the C-style indexing
       * This is an expensive operation. It is making a copy
       * of the subsequence.
       * 
       * name of new sequence will be 
       *    parent.name + "sub" + begin_end 
       *    where begin and end are 1-based index.
       * title of new sequence will be parent.title + "b-e"
       * @param b 0-based index in the sequence.
       */
      bioseq subsequenceWithName(int b, int len=-1) const;
      void setName(const string &n) { name=n; }
      /**
       * extend name. No space will be added.
       */
      void appendName(const string &aux) { name += aux; }
      void setTitle(const string &t) { title=t; }
      /**
       * Add extra text to the title using space 
       * as separator.
       */
      //void appendTitle(const string &aux) { title += " " + aux; }
      /**
       * Usper specify the separator for the additional text
       * @para aux additional string to append to the end of title.
       * @param sep seprator such as ', ' or '; ', default is SPACE.
       */
      void appendTitle(const string &aux, const string &sep="");
      bool hasTitle() const { return !title.empty(); }
      /**
       * @return the sequence id or name
       */
      const string& getName() const { return name; }
      /**
       * @return the sequence id or name
       * Alias for getName
       */
      const string& getId() const { return name; }
      const string& getTitle() const { return title; }
      /** 
       * Simply return the const reference of the underlying string.
       * Note: you cannot modify the underlying sequence.
       * The toString() function return a copy.
       * @return the sequence as string reference
       *   You cannot modify the string.
       */
      const string& getSequence() const { return seq; }

      /**
       * return the C-style index of the subseq location.
       *    -1 for not found
       */
      int locateSubsequence(const string &subseq) const;
      /**
       * Will only look for patter between b and e inclusive [b, e] c_style
       * index. Even the pattern is just on e, it will return.
       */
      int locateSubsequence(const string &subseq, const size_t b, const size_t e) const;
      int locateSubsequence(const string &subseq, const size_t b) const;

      /** the underlying sequence has been changed!
       * the title has also been updated
       */
      void rvc();
      /** 
       * retern all letter frequency regardless of case
       * sum of all frquency is 1.
       */
      map<char, double> getFrequency() const;

      /** 
       * compute the entropy of this sequence
       * for both single- and double-letter 
       * @return <1,2> residue entropy.
       */
      pair<double,double> computeEntropy() const {
         return entropy(seq); 
      }
      /** 
       * @param b is the 0-based index
       * @param e is the end of the index (inclusive)
       */
      pair<double,double> computeEntropy(int b, unsigned e) const {
         if (b<0 || e>seq.length()-1) { 
            cerr << "index error in computeEntropy()"; 
            throw bioseqexception("index out of range inside bioseq::computeEntropy");
            //exit(1);
         }
         return entropy(seq.substr(b,e-b+1)); }

      /** 
       * @param seq input string can be anything other than DNA, but only
       *    meaningful for DNA or RNA sequences.
       * @return the single letter and double letter
       *    entrpy of the sequence string.
       *
       * This is not a member function for convenient application on string.
       * This is mainly used for judging whether input sequence is 
       * simple repeat. For more accurate method a separate method should
       * be used.
       */
      static pair<double,double> entropy(const string &seq);
      /**
       * @return the length of the underlying sequence.
       */
      size_t length() const { return seq.length(); }

      /**
       * Derived classes should either overwrite it or 
       * define different functions.
       *
       * This is a more generic methods that should work
       * for any string of characters
       * The algorithm is aachar2num(Residue) essential toupper(R) - A so A or a will be 0
       * B or b will be 1, ....
       *
       * Protein getcode also use aa2num so there is no need for protein
       * subclass to have this function.
       *
       * and Nucleic acids should use hashbase.
       * These different code will be used by different score matrices 
       * for indexing of the matrix.
       * This is important for sequence alignment.
       *
       * This operation also set the local storge.
       * The code filed is not filled by default for performance
       * reasons.
       *
       * The default implementaton is the amino acid version.
       * So the Protein class does not overrite this method.
       */
      virtual const int* getcode() const;

      /** return the ith residue. Residue at index i
       * 0-based index.
       */
      char operator[](const int i) const { return seq[i]; }
      /**
       * @return a reference to the ith residue
       */
      char& operator[](const int i) { return seq[i]; }
      /**
       * Assignment operator.
       * Will not copy code incase you are construcing a different
       * object.
       */
      bioseq& operator=(const bioseq &s);
      /**
       * Move assignment operator
       */
      bioseq& operator=(bioseq &&s);
      bioseq& operator+=(const bioseq &s) { seq += s.seq; return *this; }
      /**
       * combine this and another sequence
       */
      bioseq operator+(const bioseq &s) const;

      bioseq& operator=(const string &str);
      /**
       * The underlying sequence is the same, name and title are ignored in the
       * comparison
       */
      bool operator==(const bioseq &s) const { return seq == s.toString(); }
      bool operator<(const bioseq &s) const { return seq < s.toString(); }
      /** test whether the sequence is empty or not */
      bool empty() { if (seq.empty()) return true; else return false; }
      void randomize() { 
         random_shuffle(seq.begin(), seq.end()); 
         delete[] code;
         code=nullptr; 
      }
      /** 
       * convert the underlying sequence into all lower case letters */
      void toLowerCase();
      /** 
       * convert the underlying sequence into all upper case letters */
      void toUpperCase();
      virtual SequenceType getSequenceType() const { return GENERIC; }
      /**
       * Figure out the type of the underlying sequence.
       */
      SequenceType guessType() const;

      /// debug functions ///
      void show() const;

      static bool isBioseq(const string& str) {
         return str.find_first_of("ABCDEFGHIKLMNPQRSTUVWXYZ*abcdefghiklmnpqrstuvwxyz") != string::npos;
      }
      static bool isNotBioseq(const string& str) {
         return str.find_first_not_of("ABCDEFGHIKLMNPQRSTUVWXYZ*abcdefghiklmnpqrstuvwxyz") != string::npos; 
      }

   protected:
      /**
       * The case for seq is preserved as given.
       * For DNA and AA, it is usually written as upper case.
       * The DNA identity matrix uses Upper case. 
       */
      string seq;
      string name;
      string title;
      /** 
       * this class does nothing to the code. The derived class
       * works on this member.
       * 26 letter representation of the sequence symbols.
       * touuper(Residue One-letter code) - 'A''s ASCII code.
       * This converts the letter to 0-25 for indexing operations and look up
       * in some fast operation.
       * For DNA, I used 2-bit coding A -> 0, C -> 1, G-> 2, T -> 3, S -> 4 ...
       * See codon.h and codon.cpp for detail.
       * Protein class should not need to implement the conversion.
       *
       * This is a lazy filed that is only filled when getcode() is called.
       */
      mutable int* code;
      /**
       * exchange code. To be used by revcomp series
       */
      void swapCode(const bioseq& s) { swap(code, s.code); }
};


/** all protein single letter code use upper case
 * The most recent version of amino acid code is stored in this class.
 */
class Protein : public bioseq {
   public:
      Protein() : bioseq() { codefreq[0]=-1; }
      /**
       * Copy constructor.
       */
      Protein(const Protein &s) : bioseq(s), codefreq()
      { 
         for (size_t i=0; i<numcodes; i++) codefreq[i]=s.codefreq[i]; 
      }
      /**
       * Move constructor.
       */
      Protein(Protein &&p) : bioseq(std::move(p)),
         codefreq(std::move(p.codefreq))
      {
      }

      Protein(const string &str) : bioseq(str) { codefreq[0]=-1; }
      Protein(const bioseq &s) : bioseq(s)  { codefreq[0]=-1; }
      Protein(const string &n, const string &s) 
         : bioseq (n,s)  { codefreq[0]=-1; }
      Protein subseq(int b, int e) const {
         return bioseq::subseq(b,e); }
      ~Protein() { /* parent destructor is enough */ }
      Protein& operator=(const Protein &p);
      Protein& operator=(Protein &&p);
      Protein& operator=(const string &str);

      static const char* threeLetterCode(const char A) 
      { return symbols[2*(toupper(A)-'A')]; }
      static const char* aminoAcidFullName(const char A) 
      { return symbols[2*(toupper(A)-'A')+1]; }

      /** relative to SwissProt statistics */
      double relativeEntropy() const;
      double relativeEntropyUniform() const;

      // protein class use the bioseq::getcode()

      /** the array is according the 26 letter codes
       * Note: all letters are used.
       */
      //const double* getCodeFrequence();
      array<double, 27>& getCodeFrequence();
      static const int numcodes=27; // including *
      static const char* oneLetterSymbols;
      /** for good programming this method is not needed
       */
      SequenceType getSequenceType() const { return PROTEINSEQ; }
      bool hasStart() const { return seq[0] == 'M'; }
      bool hasStop() const { return seq[seq.length()-1] == '*'; }
      bool hasInternalStop() const;
      int countInternalStops();

   private:
      /** stores the amino acid frequence as an array
       * A in zero, B in one ....
       * Calling computeEntropy will fill this array
       * The first element is set to -1 to indicate that
       * this array has not been set.
       * Not occuring letters have frequency of zero.
       * # 27 is the stop codon "*"
       */
      //mutable double codefreq[27];
      // the following line has type issues when used at pointer level
      mutable array<double, 27> codefreq;

      /**
       * array of 2x27 three letter code, description
       * The char(index + A) is the single letter code.
       */
      static const char* symbols[];
      /** frequency taken from SwissProt */
      static const double aafreq[];
      /** assume uniform distribution of Bases,
       * and uniform distribution of codons.*/
      //static const double aafrequniform[];
};

/** this class add a few DNA specific operations
 * These operation can be applied to cDNA or genomic
 * DNA
 * If sequence is using U for T, then the 
 * caller needs to convert U->T. This reduced the
 * need for costly conversion in downstream steps.
 */
class DNA : public bioseq {
   public:
      static char toRevcompBase(const char B);

      DNA() : bioseq() { }
      DNA(const DNA &s) : bioseq(s) { }
      /**
       * Move constructor
       */
      DNA(DNA &&s) : bioseq(std::move(s))  { }
      /**
       * Construct DNA from a string. Making a copy
       * of the string.
       */
      DNA(const string &s) : bioseq(s)  { }
      /**
       * Construct from a string by moving it.
       */
      DNA(string &&s) : bioseq(std::move(s))  { }

      /**
       * @param n name of the sequence
       * @param s Nucleotide Residue string
       */
      DNA(const string &n, const string &s) 
         : bioseq(n,s)  {}
      DNA(string &&n, string &&s) 
         : bioseq(std::move(n), std::move(s))  {}
      /**
       * @param n name
       * @param s sequence string
       * @param t title or description of the sequence.
       */
      DNA(const string &n, const string &s, const string &t) 
         : bioseq(n,s,t)  { }

      /**
       * Should not copy the code because the bioseq and protein
       * use the same coding algorithm, but DNA use a different 
       * algorithm.
       */
      explicit DNA(const bioseq &s) : bioseq(s)  { }
      explicit DNA(bioseq &&s) : bioseq(std::move(s))  { }
      //~DNA() { if (code != 0) delete[] code; }
      //~DNA() { bioseq::~bioseq(); }
      /**
       * Not sure the cirtual is necessary or not
       * for derived class.
       */
      ~DNA() { /* parent destructor is enought */ }
      DNA& operator=(const DNA &s);
      DNA& operator=(DNA &&s);
      //DNA& operator+=(const DNA &s) { return bioseq::operator+=(s); }
      /** replace the underlying sequece with str
       */
      DNA& operator=(const string &str);

      /** this is one essential operationn with cDNA
       * CDS features of genomic DNA also use this feature.
       * NNN => ?, may be I should use X
       */
      Protein translate() const;
      /** 
       * use 1-based index 
       * This operation limits the operation to a subrange
       * of the DNA sequence.
       * This class use only one codon table that
       * is static.
       */
      Protein translate(int begin, int end=-1) const;
      /**
       * Nothing more than the parent method.
       * 1-based index.
       */
      DNA subseq(int b, int e) const {
         return DNA(bioseq::subseq(b,e)); }
      DNA subsequenceWithName(int b, int len=-1) const {
         return DNA(bioseq::subsequenceWithName(b, len)); }
      /**
       * Overload bioseq method to return the same type
       * @param len use -1 to signal to the end of this sequence.
       */
      DNA subsequence(int b, int len=-1) const;

      /**
       * This overwrite the bioseq getcode function.
       * This use hashbase function to convert Nucleic acids to codes.
       * integers like A=>0, C=>1, ... 
       * This method is different, use a different hashing
       * method
       */
      const int* getcode() const;

      SequenceType getSequenceType() const { return DNASEQ; }
      /** return the fraction of GC/Length of sequence
       */
      double GCContent() const;
      /** another version for compute the accumulative
       * GC content for all sequences in a file
       * Add to the counts for the four bases and N, ignoring other
       * bases
       */
      double GCContent(long int &A, long int &C, long int &G, long int &T, long int &N) const;
      /**
       * Bases has only ACGT.
       */
      bool noAmbiguous() const;
      /**
       * At least one base is not ACGT
       */
      bool ambiguous() const;
      //void revcomp() { seq=reverseComplement(seq); }
      /** 
       * use the inplace global function more efficient
       * The sequence will be modified.
       */
      virtual void revcomp();
      /**
       * Make a reverse complement copy of the original sequence.
       * The name will be changed by adding _rc to the end.
       * @return a new DNA sequence that is the reverse complement of this.
       */
      DNA revcompCopy() const;
      /** return the longest peptide sequence in the
       * + strand. set b,e to the range on DNA that
       * give rise to this translation. b,e are C style
       * index. [b,e)
       */
      Protein longestORFForward(int &b, int &e) const;
      /** more efficient version, avoiding returning
       * a long string.
       **/
      void longestORFForward(Protein &p, int &b, int &e) const;
      /** [b,e) b>e from the lower strand (-)
       */
      Protein longestORFReverse(int &b, int &e) const;
      void longestORFReverse(Protein &p, int &b, int &e) const;
      Protein longestORF(char strand, int &b, int &e) const;
      Protein longestORF(int &b, int &e) const;
      void longestORF(Protein &p, int &b, int &e) const;
      /**
       * Read one sequence from the input stream ins.
       * The title of fasta sequences should not have TAB
       * but people just use it.  This function will
       * not change the TAB to space.
       * @return true if the read obtained one record
       *      false if no more sequence left in the input stream.
       * Can be used to get all sequences from one file.
       */
      bool read(istream &ins);
      /**
       * Older version to get all sequences from one fasta file.
       */
      bool read(istream &ins, string &hd) {
         return bioseq::read(ins,hd); }
      /**
       * Read one sequence from input file
       */
      void read(const string &file) { 
         bioseq::read(file); 
      }
      /**
       * @param pos is 0-based index into the DNA sequence.
       * @return the amino acid encoded by nucleotides
       *    [idx, idx+3)
       */
      char getAminoAcid(int idx) const {
         return codontable[seq.substr(idx, 3)];
      }
      /**
       * @param idx is the start of the codon
       * @param len should be 3*n, if not will expand to 3*n
       */
      string getPeptide(int idx, int len) const {
         if (len % 3 != 0) {
            len += (3 - (len%3));
         }
         return subsequence(idx, len).translate().toString();
      }
      /**
       * @return the codon starting from idx
       */
      string getCodon(int idx) const {
         return seq.substr(idx, 3);
      }
      /**
       * @param idx is 0-based index
       * @return the coding segment starting from idx for len
       *    if len is not 3n, then will padd with extra based from 
       *    this object to the 3' side
       */
      string getCodingSegment(int idx, int len) const {
         if (len % 3 != 0) {
            len += (3 - (len%3));
         }
         return seq.substr(idx, len);
      }

      static const codon& getCodonTable() { return codontable; }
      static void setCodonTable(const int ctid) { codontable.use(ctid); }

   protected:
      /** all the DNA class use one codon table 
       * This is mostly for efficiency purposes. May have to
       * add some other function to deal with particular
       * codon tables. Even global function use this codon table
       * for translation.
       **/
      static codon codontable;
};

/** 
 * should always from 5' to 3' direction
 * The underlying U is replaced by T. Still use DNA alphabet.
 * Just add coding region.
 * This object also holds a protein object.
 */
class mRNA : public DNA {
   public: 
      mRNA() : DNA() { }
      /** @param s is the input DNA sequence for the mRNA
       */
      mRNA(const string &s);
      ~mRNA() { }

   private:
      /** start of CDS 0-based index [b,e)
       */
      int cdsb;
      /** end of CDS 0-based index 1 pass the end */
      int cdse;
      Protein prt;
};

/**
 * DNA with quality score. Similar to Fastaq, but can be use for alignment
 * algorithms.
 */
class DNAQual : public DNA {
   protected:
      /**
       * Quality score for each base.
       * The allocated memory should be at least the size
       * of the sequence.
       * Needs to add one more parameter to manage the memory.
       * If this is zero, then assuming best quality.
       * This can be used to hold fasta sequences
       * where no quality score is available.
       * TODO: convert int to unsigned char to save storeage
       */
      int *qual;
      /** This is used for efficiency.
       * the getRevcomp() method should return this object.*/
      mutable DNAQual *rc;
      /**
       * default quality score 60 for fasta files
       * that does not have quality scores.
       */
      static int defaultQ;

   public: 
      DNAQual() : DNA(), qual(0), rc(0) {}
      /**
       * Copy constructor
       * Will not copy the rc that is set to zero.
       * Only if the source object has quality info,
       * the qual is copied.
       */
      DNAQual(const DNAQual &s);
      /**
       * Move constructor
       */
      DNAQual(DNAQual &&s) : DNA(std::move(s)), qual(s.qual), rc(s.rc)  {
         s.qual=0; s.rc=0;
      }
      /** Copy constructor
       * There will be no quality info associated with the object yet.
       */
      DNAQual(const string &s) : DNA(s), qual(0), rc(0)  { }

      /**
       * @param n name of the sequence
       * @param s Nucleotide Residue string
       */
      DNAQual(const string &n, const string &s) 
         : DNA(n,s), qual(0), rc(0)  {}
      /**
       * @param t is the title for the sequence.
       * @param s DNA sequence string.
       * @param n name of DNA sequence.
       * No quality is given.
       */
      DNAQual(const string &n, const string &s, const string &t) 
         : DNA(n,s,t), qual(0), rc(0)  { }
      /**
       * Make an object from name, sequence, and score.
       * @param n name of the sequence.
       * @param s actual string sequence of the BASES.
       * @param qs quality score.
       */
      DNAQual(const string &n, const string &s, int *qs);
      DNAQual(const string &n, const string &s, const unsigned char *qs);
      /**
       * @param qs Qscore range from 0-40.
       *   The largest value is 93. The character encoding add 33.
       *   char(qs+33) is the printable char encoding.
       */
      DNAQual(const string &n, const string &s, const vector<int> &qs);

      /**
       * Should not copy the code because the bioseq and protein
       * use the same coding algorithm, but DNA use a different 
       * algorithm.
       * Construct DNAQual object from DNA object.
       */
      explicit DNAQual(const DNA &s) : DNA(s), qual(0), rc(0)  { }
      /**
       * Move constructor from DNA object
       */
      DNAQual(DNA &&s) : DNA(std::move(s)), qual(0), rc(0)  { }
      /**
       * Deallocate memory.
       * deallocate the cached rc object.
       */
      ~DNAQual();
      /**
       * Assignment operator
       */
      DNAQual& operator=(const DNAQual &dq);
      /**
       * Move assignment operator.
       */
      DNAQual& operator=(DNAQual &&dq);
      //DNA& operator+=(const DNA &s) { return bioseq::operator+=(s); }
      /** replace the underlying sequece with str
       */
      DNAQual& operator=(const string &str);
      DNAQual& operator=(const DNA &dna);

      /** this is one essential operationn with cDNA
       * CDS features of genomic DNA also use this feature.
       * @param b 1-based index
       * @param e 1-based index closed interval [b, e]
       * The names will have subrange behavior as inherited from bioseq.
       */
      DNAQual subseq(int b, int e) const;
      /**
       * @param b 0-based index
       * @param len length of the subsequence if not given assuming
       *     to the end of the sequence.
       */
      DNAQual subsequenceWithName(int b, int len=-1) const;
      //{ return bioseq::subsequenceWithName(b, len); }
      /**
       * I am using virtual functions, so this should not be needed anymore.
       */
      SequenceType getSequenceType() const { return DNASEQQUAL; }
      /** 
       * use the inplace global function more efficient
       * The sequence will be modified.
       * Will reverse the quality score too.
       */
      void revcomp();
      /**
       * Make a reverse complement copy of the original sequence.
       * The name will be changed by adding _rc to the end.
       * @return a new DNA sequence that is the reverse complement of this.
       */
      DNAQual revcompCopy() const;
      /**
       * @return a pointer to the reverse complemnet of this object.
       */
      DNAQual* getRevcomp();
      /** return the quality score
       */
      int* getQuality() const { return qual; }
      /**
       * will copy the score from input qs
       * @param qs fastq quality score integer format, must have length 
       *     longer or equal to the sequence length of this object!
       */
      void setQuality(int *qs);
      void reverseQuality();
      int qscoreAt(size_t idx) const { return qual[idx]; }
      /**
       * Overload bioseq version.
       * Print the quality score to ous, width numbers per line.
       */
      virtual ostream& print(ostream &ous, int width=80) const;
      /**
       * @return the quality at position i.
       *  if fasta sequence file, then this will return 
       *  the best score.
       */
      int Qat(int i) const { if (qual == 0) return defaultQ; else return qual[i]; }
};

/**
 * Add copy number to the DNA sequences. For convenient computation
 * in some applications.
 * The quality of the object will be the average from all sequences
 * that this sequencre represents.
 */
class DNAQualCount : public DNAQual {
   protected:
      int cnt;

   public:
      DNAQualCount() : DNAQual(), cnt(0) {}
      /**
       * Copy constructor
       */
      DNAQualCount(const DNAQualCount &s) : DNAQual(s), cnt(s.cnt)  { }
      DNAQualCount(const DNAQual &s) : DNAQual(s), cnt(1)  { }
      /**
       * Making a copy of a DNAQual object. Kind of expensive.
       * Could use move constructor.
       */
      DNAQualCount(const DNAQual &s, int c) : DNAQual(s), cnt(c)  { }
      DNAQualCount(DNAQual &&s, int c) : DNAQual(std::move(s)), cnt(c)  { }
      /**
       * Move constructor
       */
      DNAQualCount(DNAQualCount &&s) : DNAQual(std::move(s)), cnt(s.cnt)  { }
      /**
       * This is the actual useful constructor making real objects.
       * @param n name of the sequence
       * @param s actual string sequence of the BASES.
       * @param qs quality score.
       * @param c count. The number of copies of this sequence.
       * Usually you are in a rush when using this object type and not
       * caring much about the title of the sequence, so the title is 
       * usually not initlaized.
       */
      DNAQualCount(const string &n, const string &s, int *qs, int c) 
         : DNAQual(n,s,qs), cnt(c) { }
      /**
       * @param qs Qscore range from 0-40.
       *   The largest value is 93. The character encoding add 33.
       *   char(qs+33) is the printable char encoding.
       */
      DNAQualCount(const string &n, const string &s, const vector<int> &qs, int c)
         : DNAQual(n,s,qs), cnt(c) { }
      DNAQualCount(DNA &&dna, int n, const vector<int> &qs);
      DNAQualCount(const DNA &dna, int n, const vector<int> &qs);
      /**
       * construct from a DNA object.
       * qual will be set to zero.
       * cnt set to 1
       */
      DNAQualCount(const DNA &dna) : DNAQual(dna), cnt(1) { }
      DNAQualCount(DNA &&dna) : DNAQual(std::move(dna)), cnt(1) { }
      /**
       * Deallocate memory.
       * deallocate the cached rc object.
       */
      ~DNAQualCount() { }
      /**
       * Assignment operator
       */
      DNAQualCount& operator=(const DNAQualCount &dq);
      /**
       * Assigment operator from parent class.
       */
      DNAQualCount& operator=(const DNAQual&dq);
      /**
       * Move assignment operator.
       */
      DNAQualCount& operator=(DNAQualCount &&dq);

      /** this is one essential operationn with cDNA
       * CDS features of genomic DNA also use this feature.
       * @param b 1-based index
       * @param e 1-based index closed interval [b, e]
       */
      DNAQualCount subseq(int b, int e) const;
      /**
       * @param b 0-based starting index
       */
      DNAQualCount subsequenceWithName(int b, int len=-1) const;

      /** 
       * you normally don't use this one.  This is poor design.
       */
      SequenceType getSequenceType() const { return DNASEQQUAL; }
      /**
       * Make a reverse complement copy of the original sequence.
       * The name will be changed by adding _rc to the end.
       * @return a new DNA sequence that is the reverse complement of this.
       */
      DNAQualCount revcompCopy() const;
      /**
       * Return the reverse complemnet version of this object.
       */
      DNAQualCount* getRevcomp();
      //{
      //   return (DNAQualCount*)DNAQual::getRevcomp();
      //}
      int qscoreAt(size_t idx) const { return qual[idx]; }
      /**
       * Overload bioseq version
       * Should be used by the '<<' operator
       * Virtual.
       * This is for debug for human to view.
       */
      ostream& print(ostream &ous, int width=80) const;
      ostream& printFastaWithHeader(ostream &ous, unsigned int width=70) const;
      /**
       * for repeated input
       */
      friend istream& operator>>(istream& ins, DNAQualCount &sq);
      void setCount(int c) { cnt=c; }
      int getCount() const { return cnt; }
};
}

#endif
