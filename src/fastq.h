#ifndef FASTQSEQ_H
#define FASTQSEQ_H

// (c) 2002 Kemin Zhou at orapra.com

#include<string>
#include<iostream>
#include <vector>
#include <array>

using namespace std;

namespace orpara {
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
       */
      int* qual;
      /** for implementation efficiency */
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
      Fastq() : name(), desc(), seq(), qual(0), qual_len(0) {  }
      /**
       * Constructor with integer quality array
       */
      Fastq(const string &id, const string &description, const string &sequence, const int* quality)
         : name(id), desc(description), seq(sequence), 
           qual(new int[sequence.length()]),
            qual_len(sequence.length())
      { 
          for (unsigned int i = 0; i<qual_len; ++i) 
             qual[i] = *(quality + i);
      }
      /**
       * Constructor with string quality
       */
      Fastq(const string &id, const string &description, const string &sequence, const string &quality)
         : name(id), desc(description), seq(sequence), 
           qual(new int[sequence.length()]),
            qual_len(sequence.length())
      { 
          encode(quality);
      }
      /**
       * Constructor with string quality, no description field
       */
      Fastq(const string &id, const string &sequence, const string &quality)
         : name(id), desc(), seq(sequence), 
           qual(new int[sequence.length()]),
            qual_len(sequence.length())
      { 
          encode(quality);
      }

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
       */
      const string& getName() const { return name; }
      /** give a new name to this sequence
       * This is useful to replace the ugly long names from sequencers.
       * */
      void setName(const string &newname) { name=newname; }
      /**
       * @return the underlying sequence as a reference.
       */
      const string& getSequence() const { return seq; }
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
      void appendDescription(const string& extra);
      void appendTitle(const string& extra) {
         appendDescription(extra);
      }

      /**
       * Given a piece of sequence it will cut at the beginning of the site
       * if the site is found in the sequence.
       * The current object will be shortened by the right piece.
       * The mothod uses string's find method.
       * @param site DNA sequence in all CAPS, upper case.
       *      | cut here
       *  ----==Site==-----
       * @return a new Fasq object to the right of site (including site).
       */
      Fastq cutLeft(const string &site);
      Fastq cutRight(const string &site);
      /**
       * Cut sequence at pos. Producing Left and Right fragments.
       *       pos
       * ======|==========
       *       +========== Right
       * ====== Left
       * 0-based index for pos. pos is part of right fragment.
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
      bool trimLowq(const unsigned int window=5, const unsigned int cutoff=20);
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
      int* getQuality() const { return qual; }
      /**
       * @return average quality score for this sequence
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
      static int getMaxScore() { return maxScore; }
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
}
#endif
