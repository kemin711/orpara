#ifndef KMERBASE_H
#define KMERBASE_H

#include <string>
#include <cstdint>

using namespace std;

/**
 * This class should be used as base for other kmer.
 * Only minimal behavior will be defined in this class.
 * For this class K can be up to 32
 */
template <int K> class KmerBase {
   public:
      KmerBase() { }
      /**
       * Helper method to convert A,C,G,T into 0-3
       * and any other base to -1
       * This should be represented by 2-bits of information.
       */
      static unsigned int base2int(char ch) {
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
               return -1;
         }
      }
      /**
       * For building the first kmer
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
       * This method should not be used for performance.
       * Instead copy paste this piece of code for particular 
       * situations.
       */
      static uint64_t* kmerize(const string& seq) {
         uint64_t* res = new uint64_t[seq.size()-K+1];
         int i;
         uint64_t v=0;
         for (i=0; i<K; ++i) { // first kmer except the last base
            v <<= 2;
            v |= base2int(seq[i]);
         }
         for (i=0; i < seq.length()-K; ++i) { // rolling update
            v<<=2;
            v |= base2int(seq[i+K]);
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


#endif
