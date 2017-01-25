#ifndef KMERMATCH_H
#define KMERMATCH_H

// (c) 2012 Kemin Zhou at orpara.com

#include <string>
#include <vector>
#include <iostream>
#include <exception>
#include <cmath>

using namespace std;

namespace orpara {
/**
 * This will be dedicated to DNA 4 bases
 * K should either be a templae parameter or 
 * static variable to gain performance.
 * This version use static member. In most cases, all classes will use the same K
 */
class Kmer {
   private:
      static int k;
      /**
       * mask based on k. For convenience.
       */
      static unsigned int mask;

      /**
       * Not sure we should store the original sequence
       */
      string seq;
      /**
       * position from 0 to L-k+1
       * hashed kmer integer value for each positon
       * Can be eliminated in the future version.
       */
      vector<int> hashval;
      /** the locations of kmers
       * index is the kmer hash value.
       * Vector store the location of the kmer's 
       * first base's 0-based index.
       */
      vector<vector<unsigned int> > loc;
      /**
       * location of the reverse complement
       */
      vector<vector<unsigned int> > locrc;

      /**
       * DNA fragements after cuting the DNA with
       * a particular kmer
       */
      vector<vector<int> > frag;
      vector<vector<int> > fragrc;

      /**
       * Compute the kmer for all positions
       * from 0 to seq.length()-k+1
       */
      void kmer2int();
      void tabulateLoc();
      void init();

   public:

      /**
       * hashval is initialized with empty vector.
       */
      Kmer() : seq(), hashval(), loc(mask+1), locrc(mask+1), 
         frag(loc.size()), fragrc(loc.size()) {}
      /**
       * hashval is intialized with -1.
       * Pow(4,k) is the same as 2<<(2*k-1)
       */
      Kmer(const string &s)
         : seq(s), hashval(s.length()-k+1, -1), 
            loc(mask+1), locrc(mask+1), 
            frag(loc.size()), fragrc(loc.size())
      {
         init();
      }

      /**
       * fill the frag vector
       */
      void digest();

      void showLocation() const;
      void showFragment() const;

      /**
       * convert sequence to an array indexed by the
       * integer representation of the word
       * @param s input string
       * @param w word size
       */
      static vector<int> hashArray(const string &s, const int w);
      /**
       * convert base to integer
       */
      static int b2i(char ch);
      /**
       * Set the mer size (k) and build the mask
       */
      static void buildMask(int mersz);
      /**
       * helper function conversting k into mask for 
       * making the hash value of kmer.
       */
      static unsigned int computeMask(int mersz);
};

/**
 * Accumulate kmer counts from a library of sequences
 * 5 is the best for telling the direction of sequences
 * Whether the 16S is in the forward or reverse direction.
 */
class KmerCount {
   private:
      unsigned int k;
      unsigned int mask;
      /**
       * index is the hash value of kmer from [0 to mask]
       * count will the number of kmer observed from input sequences.
       */
      vector<long double> count;
      vector<long double> rccount;
      int numseq;
      long int totalkmer;
      mutable vector<double> freq, rcfreq;
      mutable bool freqdone;

   public:
      KmerCount() : k(5), mask((2<<9)-1), count(mask+1), rccount(mask+1),
         numseq(0), totalkmer(0), freq(count.size()), rcfreq(count.size()),
         freqdone(false) { }
      KmerCount(unsigned int ksz) : k(ksz), mask((2<<(2*k-1))-1),
            count(mask+1), rccount(mask+1), numseq(0), totalkmer(0),
            freq(count.size()), rcfreq(count.size()),
            freqdone(false) 
      { }
      /**
       * Accumulate kmer count
       * forgot about location
       * Kmer containing non-ACGT bases will be ignored.
       * Christopher has added this version, The 
       * c perameter does not seem to be used.
       * Look like a typo error.
       */
      void operator()(const string &seq, int c);
      /**
       * Operator taking input sequences
       * @param seq the whole sequence to be digested into kmers
       *   and accumulated into the counts.
       */
      void operator()(const string &seq);
      /**
       * for observing the current object.
       */
      friend ostream& operator<<(ostream& ous, const KmerCount &kc);
      /**
       * compute the frequency from counts.
       * You need to call this function after adding
       * more sequence to this object.
       */
      void computeFrequency() const;
      const vector<double>& getFrequency() const;
      const vector<double>& getRCFrequency() const;
      /**
       * compute the distance between this and another one
       */
      double distance(const KmerCount &kc) const;
      /**
       * To see this kmer object is from the same
       * direction or not.
       */
      bool sameDirection(const KmerCount &kc) const;
      /**
       * Alias for sameDirection()
       */
      bool sameStrand(const KmerCount &kc) const { return sameDirection(kc); }
      /**
       * Read the frequency for both directrions
       */
      void open(const string &file);
      /**
       * Save freq and rcfreq to file
       * This is only used for library files that take some time
       * to process.
       */
      void save(const string &file) const;
      unsigned int getK() const { return k; }
      /** helper function that can also be used by others
       */
      static double computeD2(const vector<double> &f1, const vector<double> &f2);
};
}
#endif
