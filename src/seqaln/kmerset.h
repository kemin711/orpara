#ifndef KMERSET_H
#define KMERSET_H

#include <unordered_set>

using namespace std;

namespace orpara {

template<int K> 
class KmerSet : KmerSet<K> {
   public:
      KmerSet() : member() { }
      /**
       * Eat one sequence and convert it to  the 
       * kmer set, add the kmers to member.
       * This object can east any number of sequences.
       */
      void eat(const string& seq);
      /**
       * Add reverse complement kmer int value
       * to member.
       */
      void addRC();
      /**
       * @return the number of kmers common in both
       *   objects.
       */
      int common(const KmerSet<K>& other) const;
      int size() const {
         return member.size();
      }

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

template<int K>
void KmerSet<K>::addRC() {
   vector<int> tmp;
   for (int h : member) {
      tmp.push_back(revcompKmerInt(h));
   }
   member.insert(tmp.begin(), tmp.end());
}

template<int K>
int KmerSet<K>::common(const KmerSet<K>& other) const {
   int res=0;
   if (size() < other.size()) {
      for (int h : member) {
         if (other.member.find(h) != other.member.end()) 
            ++res;
      }
   }
   else {
      for (int h : other.member) {
         if (member.find(h) != member.end())
            ++res;
      }
   }
   return res;
}
} // orpara namespace

#endif
