#ifndef KMERBASE_H
#define KMERBASE_H
/**
 * This class should be used as base for other kmer.
 * Only minimal behavior will be defined in this class.
 */
template <int K> class KmerBase {
   public:
      KmerBase() { }
      /**
       * Helper method to convert A,C,G,T into 0-3
       * and any other base to -1
       * This should be represented by 2-bits of information.
       */
      static int base2int(char ch) {
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

   private:
      /** mask is pow(4,K)-1 
       * This should be shared by all objects in this class.
       */
      static unsigned int mask;
};

template<int K>
unsigned int KmerBase<K>::mask = (2<<(2*K-1)) - 1;

#endif