#ifndef KMERSET_H
#define KMERSET_H

#include <unordered_set>

using namespace std;

template<int K> class KmerSet : KmerSet<K> {
   public:
      KmerSet() : member() { }
      /**
       * Eat one sequence and convert it to  the 
       * kmer set, add the kmers to member.
       * This object can east any number of sequences.
       */
      void eat(const string& seq);

   private:
      /**
       * Kmers represented as integers after converting them
       * to integers.
       */
      unordered_set<int> member;
};

template<int K>
void KmerSet<K>::eat(const string& seq) {
   // calculate kmer integer hash value for all kmers in the sequence
   size_t i;
   unsigned int v=0;
   // build the intial hash value for the first kmer.
   for (i=0; i<K; ++i) {
      v <<= 2;
      v |= base2iint(seq[i]);
   }
   for (i=0; i<seq.length()-K; ++i) {
      member.insert(v);
      v<<=2;
      v |= base2int(seq[i+K]);
      v &= mask;
   }
   member.insert(v);
}

#endif
