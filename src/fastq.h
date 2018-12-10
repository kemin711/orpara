#ifndef FASTQSEQ_H
#define FASTQSEQ_H

// (c) 2002 Kemin Zhou at orapra.com

#include<string>
#include<iostream>
#include <vector>
#include <array>
#include <cstring>

using namespace std;

namespace orpara {
/**
 * TODO: convert qual from int* to unsigned char* to save space.
 */
class Fastq {
   private:
      /**
       * is the identifier of the sequence.
       * Raw reads only have very long identifiers!
       */
      string name;
      /**
       * The rest of the name line, longer description.
       */
      string desc;
      /**
       * Fastq sequence is normally in upper case
       * This program will not do the conversion for performance.
       * We assume upper case ACGT...
       * I will add a method if this is needed.
       */
      string seq;
      /**
       * For Ion Torrent, the scheme is from 0-40
       * 0                                       40
       * !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
       * 33                                      73   
       * For efficiency this is the ASCII value
       * 33 - 126
       * ASCII value
       * Qualit score should use unsigned char* to save space
       * One int use 4 bytes, char uses only one byte
       */
      unsigned char* qual;
      //int* qual;
      /** for implementation efficiency 
       * capacity of the qual, total memory (number of bytes) allocated
       * seq.size() is the actual effective usage of the total memory.
       * */
      unsigned int qual_len; // lenght of qual buffer regardless of the seq length
      /**
       * convert cscore Score in char format into the integer format
       * This will modify the underlying object.
       * The ASCII Char => Int value conversion. 
       * For Ion Torrent Platform, this should be from 33-73 integer values.
       */
      void encode(const string &cscore);
      /**
       * Will fill the cscore string with encoded chars
       * This will read from this object.
       * Integer value from 33 to 73 converted to ASCII Char vlue
       * from ! to I
       */
      void decode(string &cscore);
      void writeQualAsString(ostream &ous) const;

      static int minScore;
      /**
       * for exploration stage. Some programmers tends to make mistakes
       * This is used to capture this mistake
       */
      static int maxScore;
      /**
       * This is the value to subtract when converting
       * AscII to interger value, and added when converting
       * interger to ASCII.
       */
      static const int conv=33;

   public:
      /**
       * Default constructor.
       * This object can be used for reading from 
       * a stream.
       */
      Fastq() : name(), desc(), seq(), qual(nullptr), qual_len(0) {  }
      /**
       * Constructor with integer quality array
       */
      Fastq(const string &id, const string &description, const string &sequence, const int* quality)
         : name(id), desc(description), seq(sequence), 
           //qual(new int[sequence.length()]),
           qual(new unsigned char[sequence.length()]),
            qual_len(sequence.length())
      { 
         const int* ptr=quality;
          for (unsigned int i = 0; i<qual_len; ++i) {
             //qual[i] = *(quality + i);
             qual[i] = static_cast<unsigned char>(*ptr++);
          }
      }
      Fastq(const string &id, const string &description, const string &sequence, const unsigned char* quality)
         : name(id), desc(description), seq(sequence), 
           qual(new unsigned char[sequence.length()]),
            qual_len(sequence.length())
      { 
         memcpy(qual, quality, length());
      }

      /**
       * Constructor using string quality
       * @param id unique identifier for the sequence.
       * @param description title for the sequence
       * @param sequence the string representation of the
       *   nucleic acids sequence.
       * @param quality string representation of the quality.
       */
      Fastq(const string &id, const string &description, const string &sequence, const string &quality)
         : name(id), desc(description), seq(sequence), 
           //qual(new int[sequence.length()]),
           qual(new unsigned char[sequence.length()]),
            qual_len(sequence.length())
      { 
          encode(quality);
      }
      /**
       * Constructor with string quality, no description field
       */
      Fastq(const string &id, const string &sequence, const string &quality)
         : name(id), desc(), seq(sequence), 
           //qual(new int[sequence.length()]),
           qual(new unsigned char[sequence.length()]),
            qual_len(sequence.length())
      { 
          encode(quality);
      }
      /**
       * @param quality dynamically allocated array of unsigned char
       *    must be the same length as sequence. This object will
       *    be responsible for managing this memory.
       */
      Fastq(string &&id, string &&sequence, unsigned char* quality)
         : name(std::move(id)), desc(), seq(std::move(sequence)), 
           qual(quality),
            qual_len(length())
      {  }
      Fastq(string &&id, string&& description, string &&sequence, unsigned char* quality)
         : name(std::move(id)), desc(std::move(description)), seq(std::move(sequence)), 
           qual(quality),
            qual_len(length())
      {  }

      /**
       * Copy constructor
       */
      Fastq(const Fastq &other);
      /**
       * Move constructor
       */
      Fastq(Fastq &&other);
      ~Fastq();
      Fastq& operator=(const Fastq &other);
      /**
       * move assignment operator. Other will become useless after this
       * operation.
       */
      Fastq& operator=(Fastq &&other);
      /**
       * Use this for reading from a file.
       */
      bool read(istream &in);
      /**
       * Write out the object in fastq format.
       * If the fastq name was not prefixed with '@' it will
       * add one. Bam readers don't seem to have this special character.
       * The end-of-line will be written to the end of record.
       */
      void write(ostream &ou) const;
      /**
       * Convinient stream output format
       */
      friend ostream& operator<<(ostream &ous, const Fastq &fq) {
         fq.write(ous); return ous; }
      /**
       * write the object in fasta format with default 70 bases per line.
       */
      void writeFasta(ostream &ous, const int width=70) const;
      /** 
       * obtain the name of the sequence as string reference
       * you cannot modify this name
       * Should not contain the '@' sign.
       * Name is the unique identifier of the sequence in the file
       *   or even may be a group of related files.
       *   This function does not know whether the name is unique or 
       *   not. It is the caller's responsibility to make sure this
       *   is true.
       * @return name of sequence.
       */
      const string& getName() const { return name; }
      /** give a new name to this sequence
       * This is useful to replace the ugly long names from sequencers.
       * */
      void setName(const string &newname) { name=newname; }
      /**
       * replace long names with a short name. This is useful
       * in reducing the memory usage for very long fastq names.
       */
      void setNumberName(const string& prefix, int id) {
         setName(prefix + to_string(id));
      }
      /**
       * @return the underlying sequence as a reference.
       */
      const string& getSequence() const { return seq; }
      void setSequence(const string& sq) { 
         seq=sq;
      }
      void setSequence(string&& sq) {
         seq=std::move(sq);
      }

      void setSequenceQuality(const string& s, const string& q) {
         seq=s; encode(q);
      }
      void setSequenceQuality(string&& s, const string& q) {
         seq=std::move(s); encode(q);
      }

      string getSubSequence(int b, int len) const { return seq.substr(b,len); }
      /**
       * @return string object starting at b (0-based index) to the end of sequence
       */
      string getSubsequence(int b, int len) const { return seq.substr(b,len); }
      string getSubSequence(int b) const { return seq.substr(b); }
      string getSubsequence(int b) const { return seq.substr(b); }
      /**
       * @return a reference to the title (description of the 
       *   fastq sequence)
       */
      const string& getDescription() const { return desc; }
      const string& getTitle() const { return desc; }
      void setTitle(const string &title) { desc=title; }
      void setDescription(const string &title) { desc=title; }
      bool hasDescription() const { return !desc.empty(); }
      /**
       * Append extra information to the title.
       * In case the title was empty, this method
       * is equvelant to set title.
       */
      void appendDescription(const string& extra, char sep=' ');
      void appendTitle(const string& extra, char sep=' ') {
         appendDescription(extra, sep);
      }

      /**
       * Given a piece of sequence it will cut at the beginning of the site
       * if the site is found in the sequence.
       * The current object will be shortened by the right piece.
       * The method uses string's find method.
       * @param site DNA sequence in all CAPS, upper case.
       *      | cut here
       *  ----==site==----->
       * @return a new Fasq object to the right of site (including site).
       *  ---- current object
       *       returned object ==Site==----->
       *  Left and right reference to the cutting point
       *  meaning cut at the left side of the site
       */
      Fastq cutLeft(const string &site);
      /**
       * Similar to cutLeft
       * Cut at the right side of the site
       *              | cut here
       *  ----==site==----->
       *  after the operation
       *  ----==site== current object
       *    returned object ---->
       * @return new object after end of site
       *   if site is found in the sequence.
       */
      Fastq cutRight(const string &site);
      /**
       * Cut sequence at pos. Producing Left and Right fragments.
       *       pos
       * ======|==========
       *       +========== Right
       * ====== Left
       * 0-based index for pos. pos is part of right fragment.
       * @return the right piece. Current object
       *    will become the left piece.
       */
      Fastq cutAt(const unsigned int pos);
      /**
       * Same as above, except put result into the tailq.
       * Pos is 0-based index.
       */
      void cutAt(const unsigned int pos, Fastq &tailq);
      /**
       * 0-based index, from idx (inclusive) and right removed.
       */
      void discardTail(const unsigned int idx);
      /**
       * From [0-idx)  removed. If idx = 3, then 0,1,2 will be gone
       */
      void discardHead(const unsigned int idx);

      /**
       * Return a sub-sequence of the total with zero-indexed coordinate.
       * e must be less than the length of the sequence. 
       * inclusive meaning of coordinates.
       */
      Fastq sub(unsigned int b, unsigned int e) const;

      unsigned int length() const { return seq.length(); }
      /**
       * Fastq has not underlying sequence
       * all other fileds should also be empty but
       * non-empty values are allowed.
       */
      bool empty() const { return seq.empty(); }
      void clear() { seq.clear(); }
      /**
       * remove low quality part of the sequence. The algorithm is simple.
       * This object stores the un-decoded values of the quality score.
       * It computes a sum of quality scores in the size of the moving window.
       * When the average score of the window is below cutoff then the sequence
       * to the right (including the window) will be discarded.
       *
       * @param cutoff is the cutoff quality score, default 20. This values is
       *     the average over the window size.
       * @param window the size of the window to use default 5.
       * @return true if sequence got trimmed.
       */
      bool trimLowq(const unsigned int window=5, const unsigned int cutoff=15);
      bool trimG();
      /**
       * First trim G, then trim low quality. Trimming is done from the 3'-end.
       * @return true if either G or low qulaity regions removed from the 3' direction
       */
      bool trimGLowq(unsigned int window=6, unsigned int cutoff=10) {
         bool a=trimG(); bool  b=trimLowq(window,cutoff);
         return a || b;
      }
      /**
       * Sequence has > 0.3 N or > 0.7 single Base
       * @param nfrac cutoff for fraction of N based to call this sequence junk.
       * @param bfrac cutoff for single base fraction to call this sequence junk.
       */
      bool plaguedBySingleBase(float nfrac=0.3, float bfrac=0.75) const;
      /**
       * Use trimGLowq() then use plaguedBySingleBase()
       * @return true if trimming happened.
       */
      bool qualityTrim(unsigned int window=6, unsigned int cutoff=10, int lencut=19);
      /**
       * @return the fraction of N bases in the sequence.
       */
      std::array<float, 5> baseFraction() const;
      /**
       * reverse complement the under object.
       * After this operation, this object will be in the opposite direction.
       */
      void revcomp();
      /**
       * @return the quality pointer.
       * Note the integer values are in 33-126 range.
       * This for efficiency. To get the subtracted value
       * you should use getQscore().  If you do average
       * etc calculation using this one, make sure you subtract
       * the value.
       */
      const unsigned char* getQuality() const { return qual; }
      void setQuality(const string& strQ) {
         encode(strQ);
      }
      /**
       * @return average quality score Q for this sequence
       */
      double getAverageQuality() const;
      /**
       * @return the range of integer value stored in quality - conv
       * This is usually from 0 to 40. for Sanger it it from 0-93
       * the underlying encoding is char 33 - char 126. This one
       * does the subtraction.
       */
      vector<int> getQscore() const;
      /**
       * shift the quality score by a number
       * Say -33 for pacbio CCS2
       */
      void shiftQuality(const int shift=-33);

      static int getMinScore() { return minScore; }
      static void setMinScore(int mins) { minScore=mins; }
      static int getMaxScore() { return maxScore; }
      static void setMaxScore(int maxs) { maxScore=maxs; }
      static int getConverter() { return conv; }
      /**
       * Q=-10*log10(p) where p is the probability of base being wrong.
       */
      static double q2p(const int qval);
      /**
       * error rate p to Phred Quality Score Q
       */
      static int p2q(const double pval);
};

} // orpara namespace end
#endif
