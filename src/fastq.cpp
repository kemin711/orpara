#include "fastq.h"
#include <exception>
#include <stdexcept>
#include <sstream>
#include "bioseq.h"
#include <utility>
#include <math.h>
#include "stddev.h"
#include <numeric>
//#include <queue>
#include <cstring>
//#define DEBUG

namespace orpara {
int Fastq::minScore=999;
int Fastq::maxScore=0;

double Fastq::q2p(const int qval) {
   return pow(10, -(double)qval/10);
}

int Fastq::p2q(const double pval) {
   return -10*log10(pval);
}

Fastq::~Fastq() {
   if (qual != 0) {
      delete[] qual;
      //qual = 0;
      //qual_len = 0;
      //seq.clear();
   }
}

Fastq::Fastq(const Fastq &other) 
   : name(other.name), desc(other.desc), seq(other.seq),
      //qual(new int[other.qual_len]),
      qual(new unsigned char[other.qual_len]),
      qual_len(other.length()) 
{
   //cerr << "Calling copy constructor: other object:\n"
   //   << other.name << " " << other.seq << " " << other.qual_len << endl;
   if (length() > 0) {
      memcpy(qual, other.qual, length());
      //for (unsigned int i=0; i<other.length(); ++i) {
      //   qual[i]=other.qual[i];
      //}
   }
}

Fastq::Fastq(Fastq &&other)
   : name(std::move(other.name)),
     desc(std::move(other.desc)),
     seq(std::move(other.seq)),
     qual(other.qual),
     qual_len(other.qual_len)
{
   other.qual = nullptr;
}

Fastq& Fastq::operator=(const Fastq &other) {
   if (this != &other) {
      name=other.name;
      desc=other.desc;
      seq = other.seq;
      if (qual == nullptr) {
         qual = new unsigned char[length()];
      }
      else if (qual_len < other.length()) {
         delete[] qual;
         qual_len = other.length();
         qual = new unsigned char[length()];
      }
      memcpy(qual, other.qual, length());
      //for (unsigned int i=0; i<length(); ++i) {
      //   qual[i]=other.qual[i];
      //}
   }
   return *this;
}

Fastq& Fastq::operator=(Fastq &&other) {
   if (this != &other) {
      name=std::move(other.name);
      desc=std::move(other.desc);
      seq = std::move(other.seq);
      qual_len = other.qual_len;
      if (qual != nullptr) {
         delete[] qual;
      }
      qual=other.qual;
      other.qual=nullptr;
   }
   return *this;
}

// convert from string to unsigned int
void Fastq::encode(const string &cscore) {
   if (qual != nullptr) {
      if (qual_len < cscore.length()) {
         delete[] qual;
         qual = new unsigned char[cscore.length()];
         qual_len = cscore.length();
      }
   }
   else {
      qual = new unsigned char[cscore.length()];
      qual_len = cscore.length();
   }
   memcpy(qual, cscore.c_str(), qual_len);
   //for (unsigned int i=0; i<cscore.length(); ++i) {
   //   qual[i] = (unsigned char)cscore[i];
   //}
}

void Fastq::decode(string &cscore) {
   cscore.resize(seq.length());
   for (unsigned int i=0; i<seq.length(); ++i) {
      cscore[i] = char(qual[i]);
   }
}

void Fastq::writeQualAsString(ostream &ous) const {
   //for (unsigned int i=0; i<seq.length(); ++i) {
   //   ous << qual[i];
   //}
   ous.write(reinterpret_cast<const char*>(qual), length());
}

// the sum of shift + qual[i] should not exceed 255
void Fastq::shiftQuality(const int shift) {
   for (size_t i=0; i<seq.length(); ++i)
      qual[i] += shift;
}

void Fastq::writeFasta(ostream &ous, const int width) const {
   ous << ">" << name;
   if (!desc.empty()) ous << " " << desc;
   ous << endl;
   unsigned int i = 0;
   while (i < seq.length()) {
      ous << seq.substr(i, width) << endl;
      i += width;
   }
}

bool Fastq::read(istream &ins) {
   string line;
   char dumy[3];
   getline(ins, line);
   int emptyCnt=0;
   while (!ins.eof() && line.empty()) {
      ++emptyCnt;
      cerr << __FILE__ << ":" << __LINE__ << ":WARN empty line inside fastq file\n";
      if (emptyCnt > 500) {
         throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) +
               ":ERROR there are 500 empty lines possible NFS misbehaving");
      }
      getline(ins, line);
   }
   if (ins.eof()) {
      return false;
   }
   string::size_type i;
   if ((i=line.find(' ')) != string::npos) {
      name = line.substr(1, i-1);
      desc = line.substr(i+1);
   }
   else {
      if (line.length() == 1) {
         throw runtime_error("name empty at line " + line);
      }
      name=line.substr(1);
      if (hasDescription()) desc.clear();
   }
   getline(ins, seq); // sequence line
   ins.read(dumy, 2); // two bytes
   if (qual == nullptr) {
      qual = new unsigned char[length()];
      qual_len=length();
   }
   else if (length() > qual_len) {
      delete[] qual;
      qual = new unsigned char[length()];
      qual_len=length();
   }
   ins.read(reinterpret_cast<char*>(qual), length()); // direc read save operation
   ins.get(); // discard \n
   return true;
}

//     | cut here
// ----==Site==
// return current object missing the left part
// current object will become the left part
Fastq Fastq::cutLeft(const string &site) {
   Fastq newfasq;
   string::size_type i;
   if ((i=seq.find(site)) != string::npos) {
      if (i == 0) { // nothing needs to be done
         return newfasq;
      }
      //cerr << "cutting " << seq << " into \n";
      newfasq.name=name + "_R";
      newfasq.seq = seq.substr(i);
      newfasq.qual_len = newfasq.length();
      newfasq.qual = new unsigned char[newfasq.length()];
      //for (unsigned int j=i; j<length(); ++j) {
      //   newfasq.qual[j-i]=qual[j];
      //}
      memcpy(newfasq.qual, qual+i, newfasq.length());
      name += "_L";
      seq.resize(i);
      //cerr << seq << " | " << newfasq.seq << endl;
   }
   return newfasq;
}
//             |-----> return this part
// ---===site==------>
// ---===site== current object
Fastq Fastq::cutRight(const string &site) {
   Fastq newfasq; // to be returned
   string::size_type i;
   if ((i=seq.find(site)) != string::npos) {
      if (i + site.length() == length()) {
         // nothing needs to be done
         return newfasq;
      }
      //cerr << "site: " << site << endl;
      //cerr << "cutting " << seq << " into \n";
      newfasq.name=name + "_R";
      newfasq.seq = seq.substr(i + site.length());
      newfasq.qual_len = newfasq.length();
      newfasq.qual = new unsigned char[newfasq.length()];
      //for (unsigned int j=i+site.length(); j<length(); ++j) {
      //   newfasq.qual[j-i-site.length()]=qual[j];
      //}
      memcpy(newfasq.qual, qual+i+site.size(), newfasq.length());
      name += "_L";
      seq.resize(i + site.length());
      //cerr << seq << " | " << newfasq.seq << endl;
   }
   return newfasq;
}

// 0   4  6
// ====|======
//     pos
//  return the right piece, this object will become the left piece
Fastq Fastq::cutAt(const unsigned int pos) {
   Fastq newfasq;
   if (pos >= 0 && pos < length()) {
      newfasq.name=name + "_R";
      newfasq.seq = seq.substr(pos);
      newfasq.qual_len = newfasq.length();
      newfasq.qual = new unsigned char[newfasq.length()];
      //for (unsigned int j=pos; j<length(); ++j) {
      //   newfasq.qual[j-pos]=qual[j];
      //}
      memcpy(newfasq.qual, qual+pos, newfasq.length());
      name += "_L";
      seq.resize(pos);
   }
   else {
      cerr << pos << " outside the range of fastq sequence with length: " 
         << length() << endl;
      //throw exception("cut index outside fastq end");
      throw range_error("cut index outside fastq end");
   }
   return newfasq;
}

void Fastq::cutAt(const unsigned int pos, Fastq &newfasq) {
   if (pos > 0 && pos < length()) {
      newfasq.name=name + "_R";
      newfasq.seq = seq.substr(pos);
      newfasq.qual_len = newfasq.length();
      newfasq.qual = new unsigned char[newfasq.length()];
      //for (unsigned int j=pos; j<length(); ++j) {
      //   newfasq.qual[j-pos]=qual[j];
      //}
      memcpy(newfasq.qual, qual+pos, newfasq.length());
      name += "_L";
      seq.resize(pos);
   }
   else {
      cerr << pos << " outside fastq sequence with length: "
         << length() << endl;
      //throw exception("cutAt index out of range");
      throw range_error("cutAt index out of range");
   }
}

void Fastq::discardTail(const unsigned int idx) {
   seq.resize(idx);
}

void Fastq::discardHead(const unsigned int idx) {
   seq = seq.substr(idx);
   for (unsigned int i=0; i+idx<seq.length(); ++i) {
      qual[i] = qual[i+idx];
   }
}

Fastq Fastq::sub(unsigned int b, unsigned int e) const {
   if (b >= 0 && e < length()) {
      ostringstream convert;
      convert << name << "_sub" << b << "_" << e;
      string subname = convert.str();

      convert.str("");
      convert.clear();
      if (!desc.empty()) {
         convert << desc << " ";
      }
      convert << "subsequence from " << b << " to " << e;
      string subdescription = convert.str();
      string subseq = seq.substr(b, e-b+1);

      return Fastq(subname, subdescription, subseq, qual+b);
   }
   else {
      cerr << "range: " << b << " - " << e << " outside fastq sequence\n";
      throw range_error("sub operation out of range");
   }
}

void Fastq::write(ostream &ous) const {
   if (name[0] != '@') ous << '@';
   ous << name;
   if (!desc.empty()) {
      if (desc[0] == '@') ous << " " << desc.substr(1);
      else ous << " " << desc;
   }
   ous << endl;
   ous << seq << endl;
   ous << "+\n";
   //writeQualAsString(ous);
   ous.write(reinterpret_cast<const char*>(qual), length());
   ous << endl;
}

bool Fastq::trimLowq(const unsigned int window, const unsigned int cutoff) {
   bool trimmed = false;
   //const int window=5;
   const int cut = cutoff*window; // single score 20
   if (seq.length() < window+2) {
      return false; // not trimming very short sequences
   }
   unsigned int i;
   int sum=0;
   for (i=0; i<window; ++i) sum += (int(qual[i]) - conv);
#ifdef DEBUG
   cerr << "base value: " << conv << endl;
   for (size_t i=0; i<length(); ++i) {
      cerr << qual[i] << '|';
   }
   cerr << endl;
   cerr << "initial sum: " << sum << " cugoff=" << cutoff
      << " window=" << window << endl;
   cerr << "Running sum: " << endl;
#endif
   for (i=0; i<seq.length()-window; ++i) {
#ifdef DEBUG
      cerr << sum << " | ";
#endif
      if (sum < cut) {
#ifdef DEBUG
         cerr << "quality at " << i << " too low\n"
            << seq.substr(i, window) << endl;
         for (int j=i; j<i+window; ++j) {
            cerr << int(qual[j])-conv << '|';
         }
         cerr << endl;
#endif
         discardTail(i);
         trimmed = true;
         break;
      }
      sum -= int(qual[i]);
      sum += int(qual[i+window]);
   }
#ifdef DEBUG
   cout << endl;
#endif
   return trimmed;
}

// if fraction(N) > nfrac and any other base fraction > 0.3
// then it has a problem.
// by composition is not very fine scale.
bool Fastq::plaguedBySingleBase(float nfrac, float bfrac) const {
   array<float, 5> baseFrac = baseFraction();
   if (baseFrac[4] > nfrac) {
      //cout << "too much N: " << fq << endl;
      for (size_t i=0; i<4; ++i) {
         if (baseFrac[i] > 0.3) {
            //cout << "too much repeats: " << fq << endl;
            return true;
         }
      }
   }
   else {
      for (size_t i=0; i<4; ++i) {
         if (baseFrac[i] > bfrac) {
            //cout << "too much repeats: " << fq << endl;
            return true;
         }
      }
   }
   return false;
}

bool Fastq::qualityTrim(unsigned int window, unsigned int cutoff, int lencut) {
   bool trimmed = trimGLowq(window, cutoff);
   if (plaguedBySingleBase()) {
      clear();
      return true;
   }
   if (trimmed) {
      if (length() < lencut) {
         clear();
      }
      return  true;
   }
   return false;
}

// trim sequence like this from the right
// TCTCTTCAGGGTCCCATGCTGGGCAGGAGGGGCCTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
// using queue drags performance
// need to replace it with better one
bool Fastq::trimG() {
   static const int window=19;
   static const int cut=17;
   int GCount=0;
   string::size_type i=length()-1, ii;
   ii=i;
   // [i, ii] is the start end of the moving window
   for (auto x=0; x<window && i > 0; ++x) {
      if (seq[i] == 'G') {
         ++GCount;
      }
      --i;
   }
   if (GCount < cut) {
      return false;
   }
   while (i > 0 && GCount >= cut) {
      if (seq[ii] == 'G') --GCount;
      if (seq[i] == 'G') ++GCount;
      --i; --ii;
   }
   if (ii < length()-1) {
      discardTail(ii);
      return true;
   }
   return false;
}

void Fastq::revcomp() {
   reverseComplementInPlace(seq);
   size_t i, L;
   L = seq.length() - 1;
   for (i=0; i < seq.size()/2; ++i) {
      swap(qual[i], qual[L-i]);
   }
   name += "rc";
   if (desc.empty())
      desc = "reverse complemented";
   else
      desc += " reverse complemented";
}

vector<int> Fastq::getQscore() const {
   vector<int> tmp(length());
   for (unsigned int i=0; i<tmp.size(); ++i) {
      tmp[i]=int(qual[i])-conv;
   }
   return tmp;
}

double Fastq::getAverageQuality() const {
   //int* qv = getQuality();
   const unsigned char* qv = getQuality();
   stddev stat;
   for (size_t i=0; i<length(); ++i) {
      stat(int(qv[i]));
   }
   return stat.getMean() - conv;
}

std::array<float, 5> Fastq::baseFraction() const {
   array<int, 5> res = {0, 0, 0, 0, 0};
   for (auto b : seq) {
      switch(b) {
         case 'A' :
         case 'a' : ++res[0];
                    break;
         case 'C' :
         case 'c' : ++res[1];
                    break;
         case 'G' :
         case 'g' : ++res[2];
                    break;
         case 'T' :
         case 't' : ++res[3];
                    break;
         case 'N' :
         case 'n' : ++res[4];
                    break;
         default: cerr << "Wrong base: " << b << endl;
                  throw runtime_error("Base not ACGTN");
      }
   }
   array<float, 5> frac;
   int total = accumulate(res.begin(), res.end(), 0);
   for (size_t i=0; i<5; ++i) {
      frac[i] = res[i]/float(total);
   }
   return frac;
}

void Fastq::appendDescription(const string& extra, char sep) {
   if (extra[0] == ' ') {
      if (desc.empty()) desc = extra.substr(1);
      else {
         desc += (string(1, sep) + extra.substr(1));
      }
   }
   else {
      if (desc.empty()) desc = extra;
      else desc += (string(1, sep) + extra);
   }
}

}
