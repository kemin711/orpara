#ifndef KMERCOUNT_H
#define KMERCOUNT_H

// (c) 2012 Kemin Zhou at orpara.com

#include "kmerhelper.h"
#include <vector>

using namespace std;

namespace orpara {
/**
 * This class is for a collection of sequences that
 * are somehow related. Such as Bacterial or fungal 
 * 16S ribosomal RNA.
 * Accumulate kmer counts from a library of sequences
 * 5 is the best for telling the direction of sequences
 * Whether the 16S is in the forward or reverse direction.
 * This is the template version of the same class in
 * kmer.h.  So far this version has not been used by any
 * code yet.
 */
template<int K=5> class KmerCount {
   private:
      //unsigned int k;
      static unsigned int mask;
      /**
       * index is the hash value of kmer from [0 to mask]
       * count will the number of kmer observed from input sequences.
       */
      vector<int> count;
      vector<int> rccount;
      int numseq;
      int totalkmer;
      mutable vector<double> freq, rcfreq;
      mutable bool freqdone;

   public:
      /**
       * Default constructor of k=5
       */
      KmerCount() : count(mask+1), rccount(mask+1),
         numseq(0), totalkmer(0), freq(count.size()), rcfreq(count.size()),
         freqdone(false) { }
      /**
       * Accumulate kmer count
       * forgot about location
       * Kmer containing non-ACGT bases will be ignored.
       */
      void operator()(const string &seq);
      /**
       * for observing the current object.
       */
      template<int K>
      friend ostream& operator<<(ostream& ous, const KmerCount<int K> &kc);
      /**
       * compute the frequency from counts.
       * You need to call this function after adding
       * more sequence to this object.
       */
      void computeFrequency() const;
      const vector<double>& getFrequency() const;
      const vector<double>& getRCFrequency() const;
      /**
       * compute the distance between this and another one
       */
      double distance(const KmerCount &kc) const;
      /**
       * To see this kmer object is from the same
       * direction or not.
       */
      bool sameDirection(const KmerCount &kc) const;
      /**
       * Alias for sameDirection()
       */
      bool sameStrand(const KmerCount &kc) const { return sameDirection(kc); }
      /**
       * Read the frequency for both directrions
       */
      void open(const string &file);
      /**
       * Save freq and rcfreq to file
       * This is only used for library files that take some time
       * to process.
       */
      void save(const string &file) const;
      unsigned int getK() const { return k; }
};

//////////// KmerCount Class Implementation /////////////

template<int K>
unsigned int kmercount<K>::mask=2<<(2*K-1)-1;

// simply accumulate kmer counts
// forgot about location
template<int K>
void KmerCount<K>::operator()(const string &seq) {
   size_t i;
   unsigned int v=0;
   for (i=0; i<K; ++i) { // first kmer
      v <<= 2;
      v |= b2i(seq[i]);
   }
   for (i=0; i<seq.length()-K; ++i) {
      if (v < count.size()) {  // ignore bad kmers
         ++count[v];
         ++rccount[mask&(~v)];
      }
      else {
         cerr << "Kmer with non ACGT base ignored: "
            << seq.substr(i, K) << endl;
      }
      v<<=2;
      v |= b2i(seq[i+k]);
      v &= mask;
   }
   if (v < count.size()) {
      ++count[v];
      ++rccount[mask&(~v)];
   }
   ++numseq;
   totalkmer += (seq.length()-K+1);
}

template<int K>
ostream& operator<<(ostream& ous, const KmerCount<K> &kc) {
   if (!kc.freqdone) kc.computeFrequency();
   ous << "k: " << K << " mask: " << kc.mask << endl
      << " count size: " << kc.count.size() << endl
      << " numseq: " << kc.numseq 
      << " total kmers: " << kc.totalkmer << endl;
   for (size_t i=0; i<kc.count.size(); ++i) {
      ous << i << '\t' << kc.freq[i] << endl;
   }
   ous << string(30, 'R') << "reverse direction\n";
   for (size_t i=0; i<kc.rccount.size(); ++i) {
      ous << i << '\t' << kc.rcfreq[i] << endl;
   }
   ous << string(60, '|') << endl;
   return ous;
}

template<int K>
void KmerCount<K>::computeFrequency() const {
   for (size_t i=0; i<count.size(); ++i) {
      freq[i]=double(count[i])/totalkmer;
      rcfreq[i] = double(rccount[i])/totalkmer;
   }
   freqdone=true;
}

template<int K>
const vector<double>& KmerCount<K>::getFrequency() const {
   if (!freqdone) {
      computeFrequency();
   }
   return freq;
}

template<int K>
const vector<double>& KmerCount<K>::getRCFrequency() const {
   if (!freqdone) {
      computeFrequency();
   }
   return rcfreq;
}

template<int K>
double KmerCount<K>::distance(const KmerCount<K> &kc) const {
   vector<double> kcf=kc.getFrequency();
   double d1=computeD2(getFrequency(), kcf);
   double d2=computeD2(getRCFrequency(), kcf);
   //cout << "forward of this to kc distance: " << d1 << endl
   //   << "reverse of this to kc distance: " << d2 << endl;
   if (d1 < d2) {
      cout << "sequence in the same direction\n";
   }
   else if (d2<d1) {
      cout << "sequence in the opposite direction\n";
   }
   else {
      cout << "difficult to decide\n";
   }
   return min(d1,d2);
}

template<int K>
bool KmerCount<K>::sameDirection(const KmerCount<K> &kc) const {
   vector<double> kcf=kc.getFrequency();
   double d1=computeD2(getFrequency(), kcf);
   double d2=computeD2(getRCFrequency(), kcf);
   if (d1 < d2) {
      //cout << "sequence in the same direction\n";
      return true;
   }
   else if (d2<d1) {
      //cout << "sequence in the opposite direction\n";
      return false;
   }
   else {
      cout << "difficult to decide direction!\n";
   }
   return true;
}

template<int K>
void KmerCount<K>::save(const string &file) const {
   ofstream ouf(file);
   if (ouf.fail()) {
      throw runtime_error("failed to open file: " + file);
   }
   ouf << "k=" << K << endl
      << "numseq=" << numseq << endl
      << "totalkmer=" << totalkmer << endl;
   copy(count.begin(), count.end(), ostream_iterator<int>(ouf, "\n"));
   ouf << "reverse count:\n";
   copy(rccount.begin(), rccount.end(), ostream_iterator<int>(ouf, "\n"));
}

template<int K>
void KmerCount<K>::open(const string &file) {
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("Failed to open kmer file: " + file);
   }
   string line;
   getline(inf, line);
   if (line.substr(0,2) != "k=") {
      throw runtime_error("not kmer dump file");
   }
   int stored_k=stoi(line.substr(2));
   if (stored_k != K) {
      throw runtime_error("stored K different from template parameter " + to_string(K));
   }
   //mask=(2<<(2*k-1))-1;
   getline(inf, line);
   assert(line.substr(0, 7) == "numseq=");
   numseq=stoi(line.substr(7));
   getline(inf, line);
   assert(line.substr(0,10) == "totalkmer=");
   totalkmer=stoi(line.substr(10));
   getline(inf, line);
   count.clear();
   rccount.clear();
   count.reserve(mask+1);
   rccount.reserve(mask+1);
   while (!inf.eof() && line != "reverse count:") {
      count.push_back(stof(line));
      getline(inf, line);
   }
   getline(inf, line);
   while (!inf.eof()) {
      rccount.push_back(stof(line));
      getline(inf, line);
   }
   computeFrequency();
   freqdone=true;
}
}
#endif
