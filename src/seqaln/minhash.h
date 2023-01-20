#ifndef MINHASH_H
#define MINHASH_H

#include <kmerbase.h>
#include <map>
#include <functional>
#include <cassert>
#include <cstring>
#include <iostream>
#include <functional>
#include <cmath>

//#define DEBUG

using namespace std;

namespace orpara {

/**
 * Sequence should be free of ambiguous bases.
 * Template parameter 
 *   K is kmer length should be < 32
 *   S sketch size should be < input string length.
 *   Usually input sequence is a whole genome sequence of
 *   a bacteria or fugus. Or human. A large input file
 *   containing large number of short reads can also be 
 *   used.
 *   BHF base to integer [0-4] hashing function default is
 *   KmerBase<K>::base2int we can use a different hashing
 *   function to get more accurate estimate of similarity.
 */
//template<int K, int S, typename BHF=typename KmerBase<K>::base2int >
//template<int K, int S, function<unsigned int(const char)>BHF=&Nucleotide::base2int >
template<int K, int S> class Sketch {
   public:
      /**
       * Default constructor of empty object.
       */
      Sketch() : b2ifunc(Nucleotide::base2int), hv(new uint64_t[S]) {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX;
      }
      /**
       * Constructor from a DNA string.
       * Use base2int() function from KmerBase class.
       * @param dna input DNA sequence must have length < S and K
       */
      Sketch(const string& dna) 
         : b2ifunc(Nucleotide::base2int), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      /**
       * Use base2intC hash function
       * A special base2int function can be provided.
       * Typically you can use alternative to base2int function
       * such as those from KmerBase, this way, you can
       * use multiple hash function to reduce the false positive 
       * rates:
       *   KmerBase::base2intC(char ch), base2intG(), and base2intT()
       *
       */
      Sketch(const string& dna, function<unsigned int(const char)> hf) 
         : b2ifunc(hf), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketch(const Sketch<K,S>& sk) 
         : b2ifunc(sk.b2ifunc), hv(new uint64_t[S]) 
      {
         memcpy(hv, sk.hv, S*sizeof(uint64_t));
      }
      Sketch(Sketch<K,S>&& sk)
         : b2ifunc(sk.b2ifunc), hv(sk.hv) 
      {
         sk.hv=nullptr;
      }
      ~Sketch() {
         if (hv != nullptr) {
            delete[] hv;
         }
      }
      Sketch& operator=(const Sketch<K,S>& sk) {
         if (this != &sk) {
            b2ifunc = sk.b2ifunc;
            memcpy(hv, sk.hv, sizeof(uint64_t)*S);
         }
         return *this;
      }
      Sketch& operator=(Sketch<K,S>&& sk) {
         if (this != &sk) {
            b2ifunc = sk.b2ifunc;
            hv=sk.hv;
            sk.hv=nullptr;
         }
         return *this;
      }
      bool operator==(const Sketch<K,S>& sk) const {
         //if (b2ifunc != sk.b2ifunc) return false;
         //Cannot compare generic functions yet
         for (int i=0; i<S; ++i) {
            if (hv[i] != sk.hv[i]) return false;
         }
         return true;
      }
      bool operator<(const Sketch<K,S>& sk) const {
         for (int i=0; i<S; ++i) {
            if (hv[i] < sk.hv[i]) return true;
            if (hv[i] > sk.hv[i]) return false;
         }
         return false;
      }
      bool operator>(const Sketch<K,S>& sk) const {
         for (int i=0; i<S; ++i) {
            if (hv[i] > sk.hv[i]) return true;
            if (hv[i] < sk.hv[i]) return false;
         }
         return false;
      }

      /**
       * Jacard Similarity (A^B)/(A+B) intersect/union
       * @return Jacard similarity
       */
      float jsim(const Sketch<K,S>& sk) const {
         int i=0,j=0;
         int same=0, total=0;
         while (i<S && j < S) {
            if (hv[i] == UINT64_MAX || sk.hv[j] == UINT64_MAX) {
               total += (max(S-i, S-j));
               i=S; j=S;
               break;
            }
            else if (hv[i] == sk.hv[j]) {
               ++same; ++total;
               ++i; ++j;
            }
            else if (hv[i] < sk.hv[j]) {
               ++total; ++i;
            }
            else {
               ++total; ++j;
            }
         }
         if (i<S) {
            total += (S-i);
         }
         if (j<S) {
            total += (S-j);
         }
         return float(same)/total;
      }
      /**
       * Debug show function
       */
      ostream& show(ostream& ous) const {
         ous << hv[0];
         for (int i=1; i<S; ++i) {
            ous << ',' << hv[i];
         }
         return ous;
      }

      void add(const string& dna) {
         //assert(dna.size() > S && dna.size() > K);
         assert(dna.size() > K);
         work(dna);
      }

   private:
      /**
       * Can work on multiple short sequences to add up to longer sketches
       */
      void work(const string& dna) {
         //assert(dna.size() > S && dna.size() > K);
         assert(dna.size() > K);
         //for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         int i;
         uint64_t v=0;
         for (i=0; i<K-1; ++i) { // first kmer except the last base
            v <<= 2;
            //v |= KmerBase<K>::base2int(dna[i]);
            v |= b2ifunc(dna[i]);
         }
         for (i=0; i < (int)dna.size()-K+1; ++i) { // rolling update
            v<<=2;
            v |= b2ifunc(dna[i+K-1]);
            v &= KmerBase<K>::mask;
            if (v < hv[0]) {
               std::memmove(hv+1, hv, (S-1)*sizeof(uint64_t));
               hv[0]=v;
            }
            else if (v < hv[S-1]) {
               int x;
               bool dupval=false;
               for (x=S-1; x > -1; --x) {
                  if (v > hv[x]) break;
                  else if (v == hv[x]) { // ignore identical kmers
                     dupval=true;
                     break;
                  }
               }
               if (!dupval) { 
                  ++x;
                  std::memmove(hv+x+1, hv+x, (S-x-1)*sizeof(uint64_t));
                  hv[x]=v;
               }
            }
         }
      }

      function<unsigned int(const char)> b2ifunc;
      /**
       * array of hash values
       */
      uint64_t* hv;
};

/**
 * The size of the sketch can be controlled at run time
 */
template<int K> class Sketchv {
   public:
      /**
       * Default constructor of empty object.
       * Sketch size is 500 by default.
       */
      Sketchv() : b2ifunc(Nucleotide::base2int), S(500), hv(new uint64_t[S]) {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX;
      }

      Sketchv(int sks) : b2ifunc(Nucleotide::base2int), S(sks), hv(new uint64_t[S]) {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX;
      }
      /**
       * Constructor from a DNA string.
       * Use base2int() function from KmerBase class.
       * @param dna input DNA sequence must have length < S and K
       */
      Sketchv(const string& dna) 
         : b2ifunc(Nucleotide::base2int), S(500), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketchv(const string& dna, int sks) 
         : b2ifunc(Nucleotide::base2int), S(sks), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      /**
       * Use base2intC hash function
       * A special base2int function can be provided.
       * Typically you can use alternative to base2int function
       * such as those from KmerBase, this way, you can
       * use multiple hash function to reduce the false positive 
       * rates:
       *   KmerBase::base2intC(char ch), base2intG(), and base2intT()
       *
       */
      Sketchv(const string& dna, function<unsigned int(const char)> hf) 
         : b2ifunc(hf), S(500), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketchv(const string& dna, function<unsigned int(const char)> hf, int sks) 
         : b2ifunc(hf), S(sks), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketchv(const Sketchv<K>& sk) 
         : b2ifunc(sk.b2ifunc), S(sk.S), hv(new uint64_t[S]) 
      {
         memcpy(hv, sk.hv, S*sizeof(uint64_t));
      }
      Sketchv(Sketchv<K>&& sk)
         : b2ifunc(sk.b2ifunc), S(sk.S), hv(sk.hv) 
      {
         sk.hv=nullptr;
      }
      ~Sketchv() {
         if (hv != nullptr) {
            delete[] hv;
         }
      }
      Sketchv& operator=(const Sketchv<K>& sk) {
         if (this != &sk) {
            b2ifunc = sk.b2ifunc;
            S = sk.S;
            memcpy(hv, sk.hv, sizeof(uint64_t)*S);
         }
         return *this;
      }
      Sketchv& operator=(Sketchv<K>&& sk) {
         if (this != &sk) {
            b2ifunc = sk.b2ifunc;
            S = sk.S;
            hv=sk.hv;
            sk.hv=nullptr;
         }
         return *this;
      }
      bool operator==(const Sketchv<K>& sk) const {
         //if (b2ifunc != sk.b2ifunc) return false;
         //Cannot compare generic functions yet
         if (S != sk.S) return false;
         for (int i=0; i<S; ++i) {
            if (hv[i] != sk.hv[i]) return false;
         }
         return true;
      }
      /**
       * Ordering is only meaninful between object of the
       * same hashing function and same Sketch size
       */
      bool operator<(const Sketchv<K>& sk) const {
         for (int i=0; i<S; ++i) {
            if (hv[i] < sk.hv[i]) return true;
            if (hv[i] > sk.hv[i]) return false;
         }
         return false;
      }
      bool operator>(const Sketchv<K>& sk) const {
         for (int i=0; i<S; ++i) {
            if (hv[i] > sk.hv[i]) return true;
            if (hv[i] < sk.hv[i]) return false;
         }
         return false;
      }

      /**
       * Jacard Similarity (A^B)/(A+B) intersect/union
       * @return Jacard similarity
       */
      float jsim(const Sketchv<K>& sk) const {
         int i=0,j=0;
         int same=0, total=0;
         while (i<S && j < S) {
            if (hv[i] == UINT64_MAX || sk.hv[j] == UINT64_MAX) {
               total += (max(S-i, S-j));
               i=S; j=S;
               break;
            }
            else if (hv[i] == sk.hv[j]) {
               ++same; ++total;
               ++i; ++j;
            }
            else if (hv[i] < sk.hv[j]) {
               ++total; ++i;
            }
            else {
               ++total; ++j;
            }
         }
         if (i<S) {
            total += (S-i);
         }
         if (j<S) {
            total += (S-j);
         }
         return float(same)/total;
      }
      double identity(const Sketchv<K>& sk) const {
         double j=jsim(sk);
         return 1 + log(2*j/(1+j))/K;
      }

      /**
       * Debug show function
       */
      ostream& show(ostream& ous) const {
         ous << hv[0];
         for (int i=1; i<S; ++i) {
            ous << ',' << hv[i];
         }
         return ous;
      }

      void add(const string& dna) {
         //assert(dna.size() > S && dna.size() > K);
         assert(dna.size() > K);
         work(dna);
      }

   private:
      /**
       * Can work on multiple short sequences to add up to longer sketches
       */
      void work(const string& dna) {
         //assert(dna.size() > S && dna.size() > K);
         assert(dna.size() > K);
         //for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         int i;
         uint64_t v=0;
         for (i=0; i<K-1; ++i) { // first kmer except the last base
            v <<= 2;
            //v |= KmerBase<K>::base2int(dna[i]);
            v |= b2ifunc(dna[i]);
         }
         for (i=0; i < (int)dna.size()-K+1; ++i) { // rolling update
            v<<=2;
            v |= b2ifunc(dna[i+K-1]);
            v &= KmerBase<K>::mask;
            if (v < hv[0]) {
               std::memmove(hv+1, hv, (S-1)*sizeof(uint64_t));
               hv[0]=v;
            }
            else if (v < hv[S-1]) {
               int x;
               bool dupval=false;
               for (x=S-1; x > -1; --x) {
                  if (v > hv[x]) break;
                  else if (v == hv[x]) { // ignore identical kmers
                     dupval=true;
                     break;
                  }
               }
               if (!dupval) { 
                  ++x;
                  std::memmove(hv+x+1, hv+x, (S-x-1)*sizeof(uint64_t));
                  hv[x]=v;
               }
            }
         }
      }

      function<unsigned int(const char)> b2ifunc;
      int S;
      /**
       * array of hash values
       */
      uint64_t* hv;
};

/**
 * The size of the sketch can be controlled at run time
 * Caller has full control of all three parameters.
 */
class Sketchkv {
   public:
      /**
       * Default constructor of empty object.
       * Sketch size is 500 by default.
       */
      Sketchkv() : b2ifunc(Nucleotide::base2int), K(13), S(500), hv(new uint64_t[S]) {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX;
      }

      Sketchkv(int kl, int sks) 
         : b2ifunc(Nucleotide::base2int), K(kl), S(sks), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX;
      }
      /**
       * Constructor from a DNA string.
       * Use base2int() function from KmerBase class.
       * @param dna input DNA sequence must have length < S and K
       */
      Sketchkv(const string& dna) 
         : b2ifunc(Nucleotide::base2int), K(13), S(500), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketchkv(const string& dna, int kl, int sks) 
         : b2ifunc(Nucleotide::base2int), K(kl), S(sks), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      /**
       * Use base2intC hash function
       * A special base2int function can be provided.
       * Typically you can use alternative to base2int function
       * such as those from KmerBase, this way, you can
       * use multiple hash function to reduce the false positive 
       * rates:
       *   KmerBase::base2intC(char ch), base2intG(), and base2intT()
       *
       */
      Sketchkv(const string& dna, function<unsigned int(const char)> hf) 
         : b2ifunc(hf), K(13), S(500), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketchkv(const string& dna, function<unsigned int(const char)> hf, int kl, int sks) 
         : b2ifunc(hf), K(kl), S(sks), hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         work(dna);
      }

      Sketchkv(const Sketchkv& sk) 
         : b2ifunc(sk.b2ifunc), K(sk.K), S(sk.S), hv(new uint64_t[S]) 
      {
         memcpy(hv, sk.hv, S*sizeof(uint64_t));
      }
      Sketchkv(Sketchkv&& sk)
         : b2ifunc(sk.b2ifunc), K(sk.K), S(sk.S), hv(sk.hv) 
      {
         sk.hv=nullptr;
      }
      ~Sketchkv() {
         if (hv != nullptr) {
            delete[] hv;
         }
      }

      Sketchkv& operator=(const Sketchkv& sk);
      Sketchkv& operator=(Sketchkv&& sk);
      bool operator==(const Sketchkv& sk) const {
         //if (b2ifunc != sk.b2ifunc) return false;
         //Cannot compare generic functions yet
         if (K != sk.K || S != sk.S) return false;
         for (int i=0; i<S; ++i) {
            if (hv[i] != sk.hv[i]) return false;
         }
         return true;
      }
      /**
       * Ordering is only meaninful between object of the
       * same hashing function and same Sketch size
       */
      bool operator<(const Sketchkv& sk) const {
         for (int i=0; i<S; ++i) {
            if (hv[i] < sk.hv[i]) return true;
            if (hv[i] > sk.hv[i]) return false;
         }
         return false;
      }
      bool operator>(const Sketchkv& sk) const {
         for (int i=0; i<S; ++i) {
            if (hv[i] > sk.hv[i]) return true;
            if (hv[i] < sk.hv[i]) return false;
         }
         return false;
      }

      /**
       * Jacard Similarity (A^B)/(A+B) intersect/union
       * @return Jacard similarity
       */
      float jsim(const Sketchkv& sk) const;
      double jdist(const Sketchkv& sk) const {
         double j=jsim(sk);
         if (j == 0) return 1;
         return -(log(2*j/(1+j))/K);
      }
      /**
       * @return estimated sequence identity from Jacard similarity
       */
      double identity(const Sketchkv& sk) const {
         return 1 - jdist(sk);
      }

      /**
       * Debug show function
       */
      ostream& show(ostream& ous) const {
         ous << hv[0];
         for (int i=1; i<S; ++i) {
            ous << ',' << hv[i];
         }
         return ous;
      }

      void add(const string& dna) {
         assert(dna.size() > static_cast<unsigned int>(K));
         work(dna);
      }
      /**
       * A hash function that convert [ACGT] to [0123]
       * in random order. There are 4 predefined functions
       * in Ncleotide class to be used.
       * Should call this function before calling the add function.
       */
      void setHash(const function<unsigned int(const char)>& hashf) {
         b2ifunc = hashf;
      }

      uint64_t makeMask() const {
         return ((~uint64_t(0))<<(64-2*K))>>(64-2*K);
      }

   private:
      /**
       * Can work on multiple short sequences to add up to longer sketches
       */
      void work(const string& dna);

      /// members ////
      function<unsigned int(const char)> b2ifunc;
      int K;
      int S;
      /**
       * array of hash values
       */
      uint64_t* hv;
};
} 
#endif
