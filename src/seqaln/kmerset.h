#ifndef KMERSET_H
#define KMERSET_H

#include <unordered_set>
#include "kmerbase.h"

using namespace std;

namespace orpara {

template<int K> 
class KmerSet : public KmerBase<K> {
   public:
      KmerSet() : member(), rcmember() {
        //cout << "KmerSet default constructor\n";
      }
      /**
       * Eat one sequence and convert it to  the 
       * kmer set, add the kmers to member.
       * This object can eat any number of sequences
       * by repeated call of this method.
       * @param seq input sequence to be converted to 
       *    the member of this object.
       */
      void eat(const string& seq);
      /**
       * Add reverse complement kmer int value
       * to member.
       */
      void addRC();
      /**
       * @param rc to use forward or reverse member of other to
       *   compute the common number reverse member of other.
       *   Default is to compute common with forward.
       * @return the number of kmers common in forward direction 
       *   which is using the member to compute.
       */
      int commonForward(const KmerSet<K>& other, const bool rc=false) const;
      /**
       * @return number of common kmers shared by reverse direction
       *   which is the rcmember
       */
      int commonReverse(const KmerSet<K>& other, const bool rc=false) const;
      /**
       * Compute for both forward and reverse to other
       * @return the larger of the shared kmers.
       */
      int common(const KmerSet<K>& other) const {
         int f = commonForward(other);
         int r = commonReverse(other);
         return  max(f,r);
      }
      /**
       * Compare both the forward and reverse kmerset of this object with
       * an input sequence's forward kmerset.
       * @return the larger of the shared kmers.
       */
      int common(const string& seq) const;

      /**
       * rcmember may not have the same size as the
       * member if use even K?
       */
      int size() const {
         return member.size();
      }

      void show(ostream& ous) const {
         ous << "member size: " << member.size() 
            << " rcmember size: " << rcmember.size() << endl;
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
   //cout << __func__ << "(" << seq << ")\n";
   // calculate kmer integer hash value for all kmers in the sequence
   size_t i;
   unsigned int v=0;
   // build the intial hash value for the first kmer.
   for (i=0; i<K; ++i) {
      v <<= 2;
      v |= KmerBase<K>::base2int(seq[i]);
   }
   for (i=0; i<seq.length()-K; ++i) {
      // skip N
      if (toupper(seq[i]) != 'N') {
         member.insert(int(v));
      }
      v<<=2;
      v |= KmerBase<K>::base2int(seq[i+K]);
      v &= KmerBase<K>::mask;
   }
   member.insert(int(v));
}

template<int K>
void KmerSet<K>::addRC() {
   for (int h : member) {
      rcmember.insert(KmerBase<K>::revcompKmerInt(h));
   }
   //member.insert(tmp.begin(), tmp.end());
}
template<int K>
int KmerSet<K>::commonForward(const KmerSet<K>& other, const bool rc) const {
   const unordered_set<int>* tmp = &other.member;
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
int KmerSet<K>::commonReverse(const KmerSet<K>& other, const bool rc) const {
   const unordered_set<int>* tmp = &other.member;
   if (rc) {
      tmp = &other.rcmember;
   }

   int res=0;
   if (rcmember.size() < tmp->size()) {
      for (int h : rcmember) {
         if (tmp->find(h) != tmp->end()) 
            ++res;
      }
   }
   else { // tmp smaller look up in rcmember
      for (int h : *tmp) {
         if (rcmember.find(h) != rcmember.end())
            ++res;
      }
   }
   return res;
}

template<int K>
int KmerSet<K>::common(const string& seq) const {
   KmerSet<K> tmp;
   tmp.eat(seq);
   return common(tmp);
}

} // orpara namespace

#endif
