#include "minhash.h"

namespace orpara {

float Sketchkv::jsim(const Sketchkv& sk) const {
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

void Sketchkv::work(const string& dna) {
   assert(dna.size() > static_cast<unsigned int>(K));
   uint64_t mask = makeMask();
   int i;
   uint64_t v=0;
   for (i=0; i<K-1; ++i) { // first kmer except the last base
      v <<= 2;
      v |= b2ifunc(dna[i]);
   }
   for (i=0; i < (int)dna.size()-K+1; ++i) { // rolling update
      v<<=2;
      v |= b2ifunc(dna[i+K-1]);
      //v &= KmerBase<K>::mask;
      v &= mask;
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

Sketchkv& Sketchkv::operator=(const Sketchkv& sk) {
   if (this != &sk) {
      b2ifunc = sk.b2ifunc;
      K = sk.K; S = sk.S;
      memcpy(hv, sk.hv, sizeof(uint64_t)*S);
   }
   return *this;
}
Sketchkv& Sketchkv::operator=(Sketchkv&& sk) {
   if (this != &sk) {
      b2ifunc = sk.b2ifunc;
      K = sk.K; S = sk.S;
      hv=sk.hv;
      sk.hv=nullptr;
   }
   return *this;
}

} // end of namespace
