#ifndef MINHASH_H
#define MINHASH_H

#include <kmerbase.h>
#include <map>
#include <functional>
#include <cassert>

//#define DEBUG

using namespace std;

namespace orpara {

/**
 * Sequence should be free of ambiguous bases.
 * Template parameter 
 *   K is kmer length should be < 32
 *   S sketch size should be < input string length
 */
template<int K, int S>
class Sketch {
   public:
      //Sketch(const string& dna) : hv{UINT64_MAX} {
      /**
       * Default constructor of empty object.
       */
      Sketch() : hv(new uint64_t[S]) {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX;
      }
      /**
       * Constructor from a DNA string.
       * Use base2int() function from KmerBase class.
       * @param dna input DNA sequence must have length < S and K
       */
      Sketch(const string& dna) 
         : hv(new uint64_t[S]) 
      {
         assert(dna.size() > S && dna.size() > K);
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         int i;
         uint64_t v=0;
         for (i=0; i<K-1; ++i) { // first kmer except the last base
            v <<= 2;
            v |= KmerBase<K>::base2int(dna[i]);
         }
#ifdef DEBUG
         cout << "i=" << i << " preparation ready\n";
         map<string, int> kmercnt;
         cout << __LINE__ << ":DEBUG last kmer idx=" << dna.size()-K+1 << " dan len=" << dna.size() << endl;
         string tmpk;
#endif
         for (i=0; i < (int)dna.size()-K+1; ++i) { // rolling update
            //tmpk=dna.substr(i,K);
            //++kmercnt[dna.substr(i,K)];
            //++kmercnt[tmpk];
            v<<=2;
            v |= KmerBase<K>::base2int(dna[i+K-1]);
            v &= KmerBase<K>::mask;
            //cout << tmpk << " i=" << i << " v=" << v << endl;
            if (v < hv[0]) {
               memmove(hv+1, hv, (S-1)*sizeof(uint64_t));
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
                  memmove(hv+x+1, hv+x, (S-x-1)*sizeof(uint64_t));
                  hv[x]=v;
               }
            }
         }
         //cout << "last i=" << i << " K=" << K << endl;
         //cout << "last i=" << i << " " << kmercnt.size() << " Kmers for this seq:\n";
         //for (auto& kc : kmercnt) {
         //   cout << kc.first << " -> " << kc.second << endl;
         //}
         //cout << string(50, '-') << endl;
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
      Sketch(const string& dna, function<unsigned int(char)> hf) 
         : hv(new uint64_t[S]) 
      {
         for (int x=0; x<S; ++x) hv[x]=UINT64_MAX; // TODO: make this faster?
         int i;
         uint64_t v=0;
         for (i=0; i<K-1; ++i) { // first kmer except the last base
            v <<= 2;
            //v |= KmerBase<K>::base2intC(dna[i]);
            v |= hf(dna[i]);
         }
         for (i=0; i < (int)dna.size()-K+1; ++i) { // rolling update
            v<<=2;
            //v |= KmerBase<K>::base2intC(dna[i+K-1]);
            v |= hf(dna[i]);
            v &= KmerBase<K>::mask;
            if (v < hv[0]) {
               memmove(hv+1, hv, (S-1)*sizeof(uint64_t));
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
                  memmove(hv+x+1, hv+x, (S-x-1)*sizeof(uint64_t));
                  hv[x]=v;
               }
            }
         }
      }
      Sketch(const Sketch<K,S>& sk) 
         : hv(new uint64_t[S]) 
      {
         memcpy(hv, sk.hv, S*sizeof(uint64_t));
      }
      Sketch(Sketch<K,S>&& sk)
         : hv(sk.hv) 
      {
         sk.hv=nullptr;
      }
      ~Sketch() {
         if (hv != nullptr) {
            delete[] hv;
         }
      }
      Sketch& operator=(const Sketch<K,S>& sk) const {
         if (this != &sk) {
            memcpy(hv, sk.hv, sizeof(uint64_t)*S);
         }
         return *this;
      }
      Sketch& operator=(Sketch<K,S>&& sk) const {
         if (this != &sk) {
            hv=sk.hv;
            sk.hv=nullptr;
         }
         return *this;
      }
      bool operator==(const Sketch<K,S>& sk) const {
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

   private:
      /**
       * array of hash values
       */
      //array<uint64_t,S> hv;
      uint64_t* hv;
};
} 
#endif
