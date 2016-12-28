#include "kmer.h"
#include <bitset> // for debug print out
#include <iterator>
#include <cassert>
#include <fstream>
#include <stdexcept>

namespace orpara {
int Kmer::k=4;
unsigned int Kmer::mask=255;

int Kmer::b2i(char ch) {
   switch (ch) {
      case 'A':
         return 0;
      case 'a':
         return 0;
      case 'C':
         return 1;
      case 'c':
         return 1;
      case 'G':
         return 2;
      case 'g':
         return 2;
      case 'T':
         return 3;
      case 't':
         return 3;
      default:
         //cerr << "base " << ch << " is not ACGT\n";
         //throw range_error("base " + ch + " ignored in kmer");
         return -1;
   }
}

unsigned int Kmer::computeMask(int mersz) {
   /*
   unsigned int mask=0;
   if (mersz>=16) throw out_of_range("kmer should < 16");
   for (size_t i=0; i<mersz; ++i) {
      mask <<=2;
      mask |= 3;
   }
   cout << "mask value: " << bitset<32>(mask) << endl;
   return mask;
   */
   return (2<<(2*mersz-1)) - 1;
}

/**
 * For external use
 */
vector<int> Kmer::hashArray(const string &s, const int w) {
   unsigned int v=0;
   vector<int> tmp(s.length()-w+1);
   // mask is 00011-w-111 total 64 bits
   unsigned int M=(2<<(2*w-1))-1;
   size_t i;
   for (i=0; i<w; ++i) {
      v <<=2; 
      v |= b2i(s[i]);
   }
   cout << "mask value: " << bitset<32>(M) << endl;
   for (i=0; i<s.length()-w; ++i) {
      //cout << bitset<32>(v) << endl;
      tmp[i]=v;
      v <<=2;
      v |= b2i(s[i+w]);
      v &= M;
   }
   tmp[i]=v;
   cout << "length of string: " << s.length()
      << " length of hash: " << tmp.size() << endl;
   return tmp;
}

///////// kmer class implementation ////////////

void Kmer::buildMask(int mersz) {
   // mask is 00011-w-111 total 64 bits
   k=mersz;
   /*
   mask=0;
   if (k>=16) throw out_of_range("kmer should < 16");
   for (size_t i=0; i<k; ++i) {
      mask <<=2;
      mask |= 3;
   }
   */
   mask=(2<<(2*k-1)) - 1;
   cout << "mask value: " << bitset<32>(mask) << endl;
}

void Kmer::kmer2int() {
   // calculate kmer integer hash value for all kmers in the sequence
   size_t i;
   unsigned int v=0;
   // build the intial hash value for the first kmer.
   for (i=0; i<k; ++i) {
      v <<= 2;
      v |= b2i(seq[i]);
   }
   for (i=0; i<seq.length()-k; ++i) {
      hashval[i]=v;
      v<<=2;
      v |= b2i(seq[i+k]);
      v &= mask;
   }
   hashval[i]=v;
}

void Kmer::tabulateLoc() {
   for (size_t i=0; i<hashval.size(); ++i) {
      if (hashval[i] != -1) {
         loc[hashval[i]].push_back(i);
      }
   }
   // convert to rc position
   //cerr << "loc size: " << loc.size() << " locrc size: "
   //   << locrc.size() << endl;
   for (unsigned int i=mask; int(i) > -1; --i) {
      if (!loc[i].empty()) {
         for (unsigned int j=loc[i].size()-1; int(j)> -1; --j) {
            locrc[mask&(~i)].push_back(seq.length()-k+1-loc[i][j]);
         }
      }
   }
}

void Kmer::init() {
   kmer2int();
   tabulateLoc();
   digest();
}

void Kmer::showLocation() const {
   for (size_t i=0; i<loc.size(); ++i) {
      cout << i << ":\n    ";
      for (size_t j=0; j<loc[i].size(); ++j) {
         cout << loc[i][j] << " ";
      }
      cout << endl;
   }
   cout << "location on the reverse strand\n";
   for (size_t i=0; i<locrc.size(); ++i) {
      cout << i << ":\n    ";
      for (size_t j=0; j<locrc[i].size(); ++j) {
         cout << locrc[i][j] << " ";
      }
      cout << endl;
   }
}

void Kmer::digest() {
   for (size_t i=0; i<loc.size(); ++i) {
      int previous=0;
      frag[i].clear();
      for (size_t j=0; j<loc[i].size(); ++j) {
         frag[i].push_back(int(loc[i][j])-previous);
         previous=loc[i][j];
      }
      frag[i].push_back(int(seq.length())-previous);
   }
   for (size_t i=0; i<locrc.size(); ++i) {
      int previous=0;
      fragrc[i].clear();
      for (size_t j=0; j<locrc[i].size(); ++j) {
         fragrc[i].push_back(int(locrc[i][j])-previous);
         previous=locrc[i][j];
      }
      fragrc[i].push_back(int(seq.length())-previous);
   }
}
void Kmer::showFragment() const {
   for (size_t i=0; i<frag.size(); ++i) {
      cout << i << '\t';
      copy(frag[i].begin(), frag[i].end(), ostream_iterator<int>(cout, "\t"));
      cout << endl;
   }
   cout << "fragments on the reverse complement\n";
   for (size_t i=0; i<fragrc.size(); ++i) {
      cout << i << '\t';
      copy(fragrc[i].begin(), fragrc[i].end(), ostream_iterator<int>(cout, "\t"));
      cout << endl;
   }
}

//////////// KmerCount Class /////////////

// simply accumulate kmer counts
// forgot about location
void KmerCount::operator()(const string &seq, int c) {
   //cerr << "string:\n" << seq << endl;
   //cerr << "length: " << seq.length() << endl;
   size_t i;
   unsigned int v=0;
    for (i=0; i<k; ++i) { // first kmer
      v <<= 2;
      v |= Kmer::b2i(seq[i]);
   }
   for (i=0; i<seq.length()-k; ++i) {
      if (v < count.size()) {  // ignore bad kmers

         ++count[v];
         ++rccount[(int)(mask&(~v))];
         ++totalkmer;
      
      }
      else {
       //  cerr << "at " << i << " Kmer with non ACGT base ignored: "
         //   << seq.substr(i, k) << endl;
      }
      v<<=2;
      v |= Kmer::b2i(seq[i+k]);
      v &= mask;
   }
   if (v < count.size()) {
      ++count[v];
      ++rccount[mask&(~v)];
      ++totalkmer;

   }
 
   ++numseq;


}


void KmerCount::operator()(const string &seq) {
   //cerr << "string:\n" << seq << endl;
   //cerr << "length: " << seq.length() << endl;
   if (seq.size() < 5*k) {
      throw runtime_error("ERROR: input sequence: " + seq + " too short.  " + __FILE__ + ": " + to_string(__LINE__));
   }
   size_t i;
   unsigned int v=0;
   for (i=0; i<k; ++i) { // first kmer
      v <<= 2;
      v |= Kmer::b2i(seq[i]);
   }
   for (i=0; i<seq.length()-k; ++i) {
      if (v < count.size()) {  // ignore bad kmers
         ++count[v];
         ++rccount[mask&(~v)];
      }
      else {
         cerr << "at " << i << " Kmer with non ACGT base ignored: "
            << seq.substr(i, k) << endl;
      }
      v<<=2;
      v |= Kmer::b2i(seq[i+k]);
      v &= mask;
   }
   if (v < count.size()) {
      ++count[v];
      ++rccount[mask&(~v)];
   }
   ++numseq;
   totalkmer += (seq.length()-k+1);
}

ostream& operator<<(ostream& ous, const KmerCount &kc) {
   if (!kc.freqdone) kc.computeFrequency();
   ous << "k: " << kc.k << " mask: " << kc.mask << endl
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

void KmerCount::computeFrequency() const {
  // totalkmer=1;
   for (size_t i=0; i<count.size(); ++i) {
      freq[i]= (count[i]);
      rcfreq[i] = (rccount[i]);

     
   }

   freqdone=true;
}

const vector<double>& KmerCount::getFrequency() const {
   if (!freqdone) {
      computeFrequency();
   }
   return freq;
}

const vector<double>& KmerCount::getRCFrequency() const {
   if (!freqdone) {
      computeFrequency();
   }
   return rcfreq;
}

double KmerCount::computeD2(const vector<double> &f1, const vector<double> &f2) {
   double d=0;
   for (size_t i=0; i<f1.size(); ++i) {
      d += pow(f1[i]-f2[i], 2);
   }
   return sqrt(d);
}

double KmerCount::distance(const KmerCount &kc) const {
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

bool KmerCount::sameDirection(const KmerCount &kc) const {
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

void KmerCount::save(const string &file) const {
   ofstream ouf(file);
   if (ouf.fail()) {
      throw runtime_error("failed to open file: " + file);
   }
   ouf << "k=" << k << endl
      << "numseq=" << numseq << endl
      << "totalkmer=" << totalkmer << endl;
   copy(count.begin(), count.end(), ostream_iterator<int>(ouf, "\n"));
   ouf << "reverse count:\n";
   copy(rccount.begin(), rccount.end(), ostream_iterator<int>(ouf, "\n"));
}

void KmerCount::open(const string &file) {
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("Failed to open kmer file: " + file);
   }
   string line;
   getline(inf, line);
   if (line.substr(0,2) != "k=") {
      throw runtime_error("not kmer dump file");
   }
   k=stoi(line.substr(2));
   mask=(2<<(2*k-1))-1;
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

