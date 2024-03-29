#ifndef KMERBASE_H
#define KMERBASE_H

#include <string>
#include <cstdint>
#include <limits>
#include <iostream>

using namespace std;

namespace orpara {

/**
  * Basic operation with nucleotide bases.
  */
class Nucleotide {
   public:
      /**
       * Helper method to convert A,C,G,T into 0-3
       * and any other base to 0 (A).
       * This should be represented by 2-bits of information.
       * For speed, we remove the exception checking.
       * Caller must make sure only 4 canoical bases in the input.
       */
      static unsigned int base2int(const char ch) {
         switch (ch) {
            case 'A': 
            case 'a':
               return 0;
            case 'C':
            case 'c':
               return 1;
            case 'G':
            case 'g':
               return 2;
            case 'T':
            case 't':
               return 3;
            default:
               //cerr << __FILE__ << ":" << __LINE__ << ":WARN non ACGT base: " << ch << endl;
               //if (!isalpha(ch)) {
               //   throw logic_error("invisible Base error");
               //}
               return 0;
               //return numeric_limits<unsigned int>::max();
         }
      }
      static const char int2base[4];
      /**
       * C hash to 0, A 1, G hash to 2, T 3
       */
      static unsigned int base2intC(const char ch) {
         switch (ch) {
            case 'A': 
            case 'a':
               return 1;
            case 'C':
            case 'c':
               return 0;
            case 'G':
            case 'g':
               return 2;
            case 'T':
            case 't':
               return 3;
            default:
               //cerr << "warn non ACGT base: " << ch << endl;
               return 0;
               //return numeric_limits<unsigned int>::max();
         }
      }
      static const char int2baseC[4];
      /**
       * C hash to 2, G hash to 0, A 1, T 3
       */
      static unsigned int base2intG(const char ch) {
         switch (ch) {
            case 'A': 
            case 'a':
               return 1;
            case 'C':
            case 'c':
               return 2;
            case 'G':
            case 'g':
               return 0;
            case 'T':
            case 't':
               return 3;
            default:
               //cerr << "warn non ACGT base: " << ch << endl;
               return 0;
               //return numeric_limits<unsigned int>::max();
         }
      }
      static const char int2baseG[4];
      /**
       * C hash to 3, G hash to 1, A 2, T 0
       */
      static unsigned int base2intT(const char ch) {
         switch (ch) {
            case 'A': 
            case 'a':
               return 1;
            case 'C':
            case 'c':
               return 2;
            case 'G':
            case 'g':
               return 3;
            case 'T':
            case 't':
               return 0;
            default:
               //cerr << "warn non ACGT base: " << ch << endl;
               return 0;
               //return numeric_limits<unsigned int>::max();
         }
      }
      static const char int2baseT[4];
};

/**
 * This class should be used as base for other kmer.
 * Only minimal behavior will be defined in this class.
 * For this class K can be up to 32. If K < 16, then
 * can cast integer to uint32_t.
 */
template <int K> class KmerBase : public Nucleotide {
   public:
      KmerBase() { }
      /**
       * Given a nucleotide string of length K, this
       * function will convert it to a 64 bit integer.
       * For building the first kmer and can also be used
       * for hashing.
       * @param kmerseq is a ACGT string of length K
       * @return has value of input K-mer using base2int() function.
       */
      static uint64_t oneMer2int(const string& kmerseq) {
           // calculate kmer integer hash value for all kmers in the sequence
         int i;
         uint64_t v=0;
         // build the intial hash value for the first kmer.
         for (i=0; i<K; ++i) {
            v <<= 2;
            //v |= b2i(seq[i]);
            v |= base2int(kmerseq[i]);
         }
         return v;
      }

      /**
       * K mer to string.
       * @param merval the integer value of a kmer.
       */
      static string oneMer2String(uint64_t merval) {
         string tmp(K, 'A'); // fill with A
         int i=K-1;
         while (i > -1) {
            unsigned int code = (merval & 3U);
            if (code == 1U) {
               tmp[i]='C';
            }
            else if (code == 2U) {
               tmp[i]='G';
            }
            else if (code == 3U) {
               tmp[i]='T';
            }
            //else {
            //   throw invalid_argument("bad character code: " + to_string(code));
            //}
            merval >>= 2;
            --i;
         }
         return tmp;
      }

      /**
       * K-1 mer to string
       */
      static string oneM1mer2String(uint64_t merval) {
         string tmp(K-1, 'A');
         int i=K-2;
         while (i > -1) {
            unsigned int code = (merval & 3U);
            if (code == 1U) {
               tmp[i]='C';
            }
            else if (code == 2U) {
               tmp[i]='G';
            }
            else if (code == 3U) {
               tmp[i]='T';
            }
            //else {
             //  throw invalid_argument("bad character code: " + to_string(code));
            //}
            merval >>= 2;
            --i;
         }
         return tmp;
      }

      /**
       * This method should not be used for performance.
       * Instead copy paste this piece of code for particular 
       * situations.
       */
      static uint64_t* kmerize(const string& seq) {
         uint64_t* res = new uint64_t[seq.size()-K+1];
         unsigned int i;
         uint64_t v=0;
         for (i=0; i<K-1; ++i) { // first kmer except the last base
            v <<= 2;
            v |= base2int(seq[i]);
         }
         for (i=0; i < seq.length()-K+1; ++i) { // rolling update
            v<<=2;
            v |= base2int(seq[i+K-1]);
            v &= mask;
            res[i]=v;
         }
         return res;
      }

      /**
       * @return the int value of the reverse complement of the kmer
       * whose int value is hv.
       */
      /*
      static unsigned int revcompKmerInt(unsigned int hv) {
         unsigned int hvc = (~hv);
         unsigned int res=(3&hvc);
         for (size_t i=0; i<K-1; ++i) {
            res <<= 2;
            hvc >>= 2;
            res |= (3&hvc);
         }
         return res&mask;
      }
      */
      static uint64_t revcompKmerInt(uint64_t hv) {
         uint64_t hvc = (~hv);
         uint64_t res=(3&hvc);
         for (size_t i=0; i<K-1; ++i) {
            res <<= 2;
            hvc >>= 2;
            res |= (3&hvc);
         }
         return res&mask;
      }

      static uint64_t getMask() { return mask; }
      static uint64_t getLeftMask() { return maskL; }
      static uint64_t getRightMask() { return maskR; }

   public:
      /** 
       * mask is pow(4,K)-1 
       * 1000...0 63 0's
       * This should be shared by all objects in this class.
       * mask is 2*K 1's to the right
       * K=3  --6---
       * 00000111111
       */
      static const uint64_t mask;
      /**
       * K-1 left word from kmer
       * For K=3
       * 111100
       */
      static const uint64_t maskL;
      /**
       * K-1 right word from kmer
       * For K=3
       * 001111
       */
      static const uint64_t maskR;
};

template<int K>
const uint64_t KmerBase<K>::mask = ((~uint64_t(0))<<(64-2*K))>>(64-2*K);
template<int K> const uint64_t KmerBase<K>::maskL=KmerBase<K>::mask << 2;
template<int K> const uint64_t KmerBase<K>::maskR=KmerBase<K>::mask >> 2;
//template<int K> const char KmerBase<K>::int2base[4]={'A', 'C', 'G', 'T'};
//template<int K> const char KmerBase<K>::int2baseC[4]={'C', 'A', 'G', 'T'};
//template<int K> const char KmerBase<K>::int2baseG[4]={'G', 'A', 'C', 'T'};
//template<int K> const char KmerBase<K>::int2baseT[4]={'T', 'A', 'C', 'G'};
} // end of namespace

#endif
