#ifndef KMERSET_H
#define KMERSET_H

#include <unordered_set>

using namespace std;

namespace orpara {

template<int K> 
class KmerSet : KmerSet<K> {
   public:
      KmerSet() : member(), rcmember() { }
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
      /**
       * rcmember may not have the same size as the
       * member if use even K?
       */
      int size() const {
         return member.size();
      }

   private:
      /**
       * Kmers represented as integers after converting them
       * to integers.
       */
      unordered_set<int> member;
      unordered_set<int> rcmember;
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
   for (int h : member) {
      rcmember.push_back(revcompKmerInt(h));
   }
   //member.insert(tmp.begin(), tmp.end());
}
template<int K>
int KmerSet<K>::commonForward(const KmerSet<K>& other, const bool rc=false) const {
   const unordered_set* tmp = &other.member;
   if (rc) {
      tmp = &other.rcmember;
   }

   int res=0;
   if (member.size() < tmp->size()) {
      for (int h : member) {
         if (tmp->find(h) != tmp->end()) 
            ++res;
      }
   }
   else {
      for (int h : *tmp) {
         if (member.find(h) != member.end())
            ++res;
      }
   }
   return res;
}

template<int K>
int KmerSet<K>::commonReverse(const KmerSet<K>& other, const bool rc=false) const {
   const unordered_set* tmp = &other.member;
   if (rc) {
      tmp = &other.rcmember;
   }

   int res=0;
   if (rcmember.size() < tmp->size()) {
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
