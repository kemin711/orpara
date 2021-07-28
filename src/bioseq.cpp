#include "bioseq.h"
#include <cctype>
//#include <map>
#include <fstream>
#include <math.h>
#include "strformat.h"

// for performance we can turn off the verification of 
// each residue, this can be costly in large inputs
//#define VERIFY_RESIDUE
using namespace std;

namespace orpara {
//using namespace KZUtility;
//////// helper function ////////////////////////////////////////
void translate(string &pep, const string &seq, unsigned int begin, unsigned int end) {
   try {
      if (!pep.empty()) pep.clear();
      codon ct=DNA::getCodonTable();
      if (begin < 1) {
         cerr << "Starting position " << begin 
            << " less than 1, this function translate() use 1-based index\n";
         exit(1);
      }
      unsigned int i=begin-1;
      if (end == 0 || end > seq.length()) 
         end = seq.length();
      while (i< end) {
         pep += ct[seq.substr(i,3)];
         i += 3;
      }
   }
   catch (exception &err) {
      cerr << "Faile to translate RNA " << seq 
         << "\nfrom " << begin << " to " << end 
         << " into peptide seq\n";
      cerr << err.what() << endl;
      exit(1);
   }
}
/*
void translate(string &pep, const string &seq, int begin, int end) {
   translate(pep, seq, (unsigned)begin, (unsigned)end); 
}
*/

int countInternalStops(const string &seq) {
   unsigned int ns=0;
   for (unsigned int i=0; i<seq.length()-1; i++) {
      if (seq[i] == '*') ++ns;
   }
   return ns;
}

char aanum2char(int c) {
   if (c == 26) return '*';
   else {
      c += 'A';
      return c;
   }
}

/** this translated full ORF for prefix and suffix ORF
  * you need special treatment
  * prefix is ok, will use the start to encode frame
  * information; if start < 3 then it contain frame
  * info.
  */
Range P2Rindex(const Range &ofb, const int frame) {
   return Range(frame+3*ofb.begin(), frame+3*ofb.end()+2);
}
Interval P2Rindex(const Interval &ofb, const int frame) {
   return Interval(frame+3*ofb.begin(), frame+3*ofb.end()+2);
}
int P2Rindex(const int pos, const int frame)
{ return 3*pos+frame; }

/** helper function for findAllORFIndex() */

bool overlayRVR(const Range &r, const vector<Range> &vr) {
   for (unsigned int i=0; i<vr.size(); ++i) {
      if (r.overlay(vr[i]) > -2) return true;
   }
   return false;
}
bool farRVR(const Range &r, const vector<Range> &vr, 
      const int hhc, const int htc, const int ttc) {
   int topo, dist;
   for (unsigned int i=0; i<vr.size(); ++i) {
      topo=vr[i].topology(r);
      dist=vr[i].distance(r);
      if ((topo == 1 && dist  < htc)
            || (topo == 2 && dist < hhc)
            || (topo == 3 && dist < ttc))
      {
         return false;
      }
   }
   return true;
}

/** @return non-overlapping ORFs that are supposed to be more 
  * biological. some of the ORF might not be the longest,
  * but the overall Coding is optimal. The algorithm is simply
  * a packing algorithm by selecting the largest first,
  * then select the next non-overlapping ORF.
  * You need to do your own filtering.
  * The result is sorted from 5' to 3' direction if rna.
  * This sorted order is essential for ESTAssembly::breakup method.
  */
vector<Range> findAllORFIndex(const string &rna, int base, 
      int peplen_cutoff, int HHcut, int HTcut, int TTcut) {
   string pep;
   int i;
   list<Range> orfbound; // tmp holder for candidate ORF
   set<Range, greaterRangeLength> allorfbound; // for sorting from large to small
   list<Range>::const_reverse_iterator rli;
   for (i=0; i<3; ++i) {
      translate(pep, rna, i+1);
      findAllPepORFIndex(orfbound, pep, peplen_cutoff);
      if (!orfbound.empty()) {
         rli=orfbound.rbegin();
         if ((unsigned)rli->end() == pep.length() - 1) { // check suffix ORF
            allorfbound.insert(Range(i+3*rli->begin(), rna.length()-1));
         }
         else allorfbound.insert(P2Rindex(*rli, i));
         ++rli;
         while (rli != orfbound.rend()) {
            allorfbound.insert(P2Rindex(*rli, i));
            ++rli;
         }
      }
   }
   string rcrna=reverseComplement(rna);
   for (i=0; i<3; ++i) {
      translate(pep, rcrna, i+1);
      findAllPepORFIndex(orfbound, pep, peplen_cutoff);
      if (!orfbound.empty()) {
         rli=orfbound.rbegin();
         if ((unsigned)rli->end() == pep.length()-1) {
            allorfbound.insert(Range(rna.length()-1-i-3*rli->begin(), 0));
         }
         else allorfbound.insert(rna.length()-1-P2Rindex(*rli, i));
         ++rli;
         while (rli != orfbound.rend()) {
            allorfbound.insert(rna.length()-1-P2Rindex(*rli, i));
            ++rli;
         }
      }
   }
   vector<Range> bestorf;
   set<Range, greaterRangeLength>::const_iterator si;
   if (allorfbound.empty()) {
      //cout << "RNA seq has no ORF:\n" << rna << endl;
      return bestorf;
   }
   // for debug, the user should inspect the result
   /*
   cout << "All orf bounds from large to small\n"
      << allorfbound.size() << " ORFs found\n";
   for (si=allorfbound.begin(); si != allorfbound.end(); ++si) {
      cout << *si << " " << si->length() << endl;
   }
   */
   si=allorfbound.begin();
   bestorf.push_back((*si)+base);
   ++si;
   while (si != allorfbound.end()) {
      //if (!overlayRVR(*si, bestorf) && farRVR(*si, bestorf, HHcut, HTcut, TTcut)) {
      if (farRVR(*si, bestorf, HHcut, HTcut, TTcut)) {
         bestorf.push_back((*si) + base);
      }
      ++si;
   }
   sort(bestorf.begin(), bestorf.end());
   // non-overlapping ORF
   /*
   cout << "\nNon-overlapping ORFs and far apart\n";
   for (i=0; i<bestorf.size(); ++i) 
      cout << bestorf[i] << " " << bestorf[i].length() << endl;
   cout << endl;
   */
   return bestorf;
}

/** find all the ORF index bounds in 0-based index
  * Number is in protein space
  * minimum AA length to register.
  * 25 aa is the minimum we are goint to register.
  * @param orfrange is the result of this operation. It will make
  *     it empty at the begining of the run if it was not.
  *     The output is in ss index. the end of the range is the 
  *     '*' symbol for complete ORF or prefix ORF. The begin
  *     of the range is the index of 'M'. 
  * @param ss is the input peptide sequence
  * @cuotff is a integer number. Peptide shorter than this are
  *     ignored in the searching phase. This parameter should
  *     not affect the real performance of this algorithm.
  *     It provided fine control. Currently I am using 25 aa.
 */
void findAllPepORFIndex(list<Range> &orfrange, const string &ss, unsigned int cutoff) {
   // Mi is the start, Si is the stop index
   if (!orfrange.empty()) orfrange.clear();
   unsigned int partialcut = max(cutoff/2, (unsigned)25);
   string::size_type Si=0;
   string::size_type Mi;
   while (Si != string::npos && (Mi = ss.find('M', Si)) != string::npos) {
      Si=ss.find('*', Mi+1);
      if (Si != string::npos) { // found stop *
         if (Si-Mi > cutoff) {
            orfrange.push_back(Range(Mi,Si));
         }
      }
      else { // no stops found
         if (ss.length()-Mi > partialcut) {
            orfrange.push_back(Range(Mi, ss.length()-1));
         }
         break;
      }
      ++Si;
   }
   // no FULL ORF found, should not have limits,
   // should always find something
   if (orfrange.empty()) {
      if (ss.length() > partialcut+1) {
         Si=ss.find('*');
         if (Si == string::npos) {
            orfrange.push_back(Range(0,ss.length()-1));
         }
         else if (Si > partialcut) {
            orfrange.push_back(Range(0,Si));
         }
      }
      return;
   }
   // if first full ORF starts after 100 there is a chance of a prefix-ORF
   if (orfrange.front().begin() > 90) {
      Si=ss.find('*');
      if (Si != string::npos && (int)Si < orfrange.front().begin()) {
         // Stop before first full ORF over 100aa,
         // and N-partial ORF longer than 35
         if (Si > partialcut && orfrange.front().begin()-Si > 90) {
            orfrange.push_front(Range(0,Si));
         }
      }
      else if (orfrange.front().begin() > orfrange.front().length()
            || (unsigned)orfrange.front().begin() > cutoff+10) {
         // --prefix orf M--ORF---* no stop between prefix and M
         // shold convert first Full ORF into a prefix ORF
         // only if prefix ORF long enough
         orfrange.front().setBegin(0);
      }
   }
}

int aachar2num(char a) {
   if (a == '*') { // stop codon special
      return 26;
   }
   else {
      int i = toupper(a) - 'A';
      return i<0? 27 : i;
   }
}

/////////////////////////////////////////////////////////////
////////////// bioseq class /////////////////////////////////

bool bioseq::read(istream &ins, string &header) {
   if (ins.eof()) return false;
   string line;
   string::size_type i;
   clear(); // make sure the sequence starts with empty
   // empty header means 1. not eof start of seqstream
   // 2. eof, end of seqstream.
   if (header.empty()) { 
      if (!ins.eof()) { // we are reading the first sequence.
         getline(ins,line);
         if (line[0] != '>' || line.length() < 2) {
            cerr << "not in proper fasta format, first line of fasta file: "
               << line << endl;
            throw bioseqexception("fasta format problem");
         }
      }
      else return false; // finished
   }
   else {
      line=header;
   }
   // now line contains the header line
   i = line.find(' ', 1);
   if (i != string::npos) {
      name=line.substr(1,i-1); // remove > in header
      title=line.substr(i+1);
   }
   else name=line;

   // read sequence part
   getline(ins,line);
   while (!ins.eof() && line[0] != '>') {
      if (line.length()>0 && 
            (line.find_first_of("ABCDEFGHIKLMNPQRSTUVWXYZ*") != string::npos || 
             line.find_first_of("abcdefghiklmnpqrstuvwxyz*") != string::npos)
            ) {
         if (!isprint(line[line.length()-1])) {
            line.erase(line.length()-1);
         }
         if (line.find_first_not_of("ABCDEFGHIKLMNPQRSTUVWXYZ*abcdefghiklmnpqrstuvwxyz") != string::npos) 
         {
            throw runtime_error("sequence " + name + " contains non-bioseq character! " + line);
         }
         seq += line;
      }
      getline(ins,line);
   }
   if (seq.empty()) {
      cerr << "Empty sequence something is wrong!\n";
      throw bioseqexception("empty sequence from sequence file");
   }
   if (ins.eof())  {
      header.clear();
   }
   else {
      header=line; 
   }
   return true;
}

// header not provided
bool bioseq::read(istream &ins) {
   if (ins.eof() || ins.peek() == EOF) return false;
   string::size_type i;
   getline(ins,title);
   if (title[0] != '>' || title.length() < 2) {
      throw bioseqexception(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR improper fasta header " + title);
   }
   if ((i=title.find(' ')) != string::npos) {
      name=title.substr(1,i-1);
      title=title.substr(i+1);
   }
   else {
      name=title.substr(1);
      title.clear();
   }
   // now read sequence part, clear code if not null
   if (code != nullptr) {
      delete[] code; code = nullptr;
   }
   getline(ins,seq); // this could be the whole aa-sequence
#ifdef VERIFY_RESIDUE
   if (isNotBioseq(seq)) {
      throw bioseqexception("biosequence contains illegal char");
   }
#endif
   if (ins.eof() || ins.peek() == '>') { // no more input, or seq shorter than one line
      return true;
   }
   string line;
   getline(ins,line);
   while (!ins.eof()) {
      if (!line.empty()) { // allowing empty lines inside sequence files
         if (!isprint(line[line.length()-1])) { // remove invisible variable
            line.resize(line.length()-1);
         }
#ifdef VERIFY_RESIDUE
         if (isNotBioseq(line)) {
            throw bioseqexception("illigal bioseq seq char in " + line);
         }
#endif
         seq += line;
      }
      if (ins.peek() == '>') return true;
      getline(ins,line);
   }
   if (seq.empty()) {
      cerr << "Empty sequence something is wrong!\n";
      throw bioseqexception("empty sequence from sequence file");
   }
   ins.peek();
   return true;
}

// this method only read one sequence from the file
// even there are multiple sequences in the file
bool bioseq::read(const string &file) {
   ifstream IN(file);
   if (IN.fail()) {
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
            + ":ERROR Failed to open input file: " + file);
   }
   // 1 read header
   getline(IN, title);
   if (title.front() != '>') {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR inproper fasta-foramted sequence\n";
      return false;
   }
   string::size_type i = title.find(' ');
   if (i == string::npos) {
      name=title.substr(1);  // first char in fasta header is '>'
      title.clear();
   }
   else {
      name=title.substr(1, i-1);
      title = title.substr(i+1);
   }
   if (!seq.empty()) seq.clear();
   if (code != nullptr) {
      delete[] code;
      code = nullptr;
   }
   string line;
   getline(IN, line);
   while (!IN.eof() && line[0] != '>') {
      if (!line.empty()) { // dealing with empty line case
         if (!isprint(line.back())) {
            line.resize(line.size()-1);
         }
         if (isNotBioseq(line)){
            cerr << line << " contain illegal sequence char\n";
            return false;
         }
         seq += line;
      }
      getline(IN, line);
   }
   if (!IN.eof()) {
      cerr << __FILE__ << ":" << __LINE__ << ":WARN more than one sequence in the file: " << file << "\n";
   }
   return true;
}

istream& operator>>(istream& ins, bioseq &sq) {
   static string line;
   getline(ins, line);
   if (line[0] != '>') {
      throw runtime_error("DNAQualCount object text representaton wrong");
   }
   size_t i=line.find(' ');
   if (i == string::npos) {
      sq.setName(line.substr(1));
      sq.title.clear();
   }
   else {
      sq.setName(line.substr(1, i-1));
      sq.setTitle(line.substr(i+1));
   }
   sq.seq.clear();
   getline(ins,line);
   while (!ins.eof() && line[0] != '>') {
      if (line.length()>0 && 
            (line.find_first_of("ABCDEFGHIKLMNPQRSTUVWXYZ*") != string::npos || 
             line.find_first_of("abcdefghiklmnpqrstuvwxyz*") != string::npos)) 
      {
         if (!isprint(line[line.length()-1])) { // remove invisible variable
            line.erase(line.length()-1);
         }
         sq.seq += line;
      }
      getline(ins,line);
   }
   if (sq.seq.empty()) {
      cerr << "Empty sequence something is wrong!\n";
      throw bioseqexception("empty sequence from sequence file");
   }
   ins.peek();
   return ins;
}

string& bioseq::importFasta(const string& fas) {
	string header;
	string::size_type i;
	i = fas.find('\n');
	header = fas.substr(0, i);
	seq.reserve(fas.size());
	seq = fas.substr(i+1);
	while ( (i=seq.find('\n')) != string::npos) {
		seq.erase(i, 1);
	}
	return seq;
}
 
string reverseComplement(const string& seq) {
	string tmp;
	tmp.reserve(seq.size()+1);
	int i = seq.size()-1;   // string::size_type non negative, trouble testing 0
	while (i >= 0) {
		//cerr << i << endl;
		if (seq[i] == 'a') tmp += 't'; // preserve case of original sequence
      else if (seq[i] == 'A') tmp += 'T'; 
		else if (seq[i] == 'c') tmp += 'g';
      else if (seq[i] == 'C') tmp += 'G';
		else if (seq[i] == 'g') tmp += 'c';
		else if (seq[i] == 'G') tmp += 'C';
		else if (seq[i] == 't') tmp += 'a';
		else if (seq[i] == 'T') tmp += 'A';
		else {
			if (isalpha(seq[i])) {
            if (islower(seq[i])) tmp += 'n';
            else //  must be upper
               tmp += 'N';
         }
			else tmp += seq[i];  // no translation
		}
		--i;
	}
	return tmp;
}

char complementChar(char ch) {
   if (ch == 'a') return 't'; // preserve case of original sequence
   if (ch == 'A') return 'T'; 
   if (ch == 'c') return 'g';
   if (ch == 'C') return 'G';
   if (ch == 'g') return 'c';
   if (ch == 'G') return 'C';
   if (ch == 't') return 'a';
   if (ch == 'T') return 'A';
   if (isalpha(ch)) {
      if (islower(ch)) return 'n';
      else return 'N';
   }
   else return ch; // no translation
}

void reverseComplementInPlace(string& seq) {
   unsigned int i;
   char a,b;
   for (i=0; i<seq.size()/2; i++) {
      a=complementChar(seq[i]);
      b=complementChar(seq[seq.size()-i-1]);
      seq[i]=b; seq[seq.size()-i-1]=a; // swap
   }
   if (seq.size() % 2 == 1) seq[i]=complementChar(seq[i]);
}

void bioseq::printFastaWithHeader(ostream &ous, unsigned int width) const {
	unsigned int i=0;
   ous << ">" << name;
   if (!title.empty()) { ous << "  " << title; }
   ous << endl;

	while (i < seq.size()) {
		ous << seq[i];
		if ((i+1)%width == 0) ous << endl;
		i++;
	}
	if (i%width != 0) ous << endl;
}
		
// no sequence name printed
void bioseq::printFasta(ostream &ous, unsigned int width) const {
	unsigned int i=0;
	while (i < seq.size()) {
		ous << seq[i];
		if ((i+1)%width == 0) ous << endl;
		i++;
	}
	if (i%width != 0) ous << endl;
}

/* helper function */
void printFasta(ostream &ous, const string& seq, unsigned int width)
{
   unsigned int i;
	for (i=0; i<seq.size(); i++) {
		ous << seq[i];
		if ((i+1)%width == 0) ous << endl;
	}
	if (i%width != 0) ous << endl;
}

// helper for << operator
ostream& bioseq::print(ostream &ous, int width) const {
   if (!name.empty()) {
      ous << ">" << name;
      if (!title.empty()) ous << " " << title;
      ous << endl;
   }
   printFasta(ous,width);
   return ous;
}

ostream& operator<<(ostream &ous, const bioseq &s) {
   return s.print(ous,80);
}

string bioseq::substr(int b, int e) const { 
   if (b < 1 || e<1) {
      cerr << "Using biologiest's coordinate system.  First index must be 1\n";
      exit(1);
   }
   if (b<=e) 
      return seq.substr(b-1,e-b+1);
   else { //if (b>e) 
      return seq.substr(e-1,b-e+1);
   }
}
string bioseq::substr(int b) const {
   if (b<1) {
      throw runtime_error("This function index start from 1, index " + to_string(b) + " out of bound!"); 
   }
   return seq.substr(b-1);
}
string bioseq::substring(const int b, const int len) const {
   if ((size_t)b >= length()) {
      cerr << "seq: " << seq << "\nb,len: " << b << ',' << len 
         << " b is greater than sequence length: " << length() << endl;
      throw bioseqexception("bioseq::substring start index: " + itos(b) 
            + " greater than seqlength");
   }
   else if (b < 0) {
      throw bioseqexception("b is negative inside bioseq::substring");
   }
   return seq.substr(b, len);
}
string bioseq::substring(const int b) const {
   if (b>=(int)length() || b < 0) {
      cerr << "seq: " << seq << endl
         << "b: " << b << endl;
      throw bioseqexception("bioseq::substring start index: " + itos(b) 
            + " out of range inside substring(b)");
   }
   return seq.substr(b);
}

// using human coordinates
bioseq bioseq::subseq(int b, int e) const { 
   if (b<1) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR: Coordinates b=" << b << " e=" << e 
         << " must use 1-based index\n";
      throw bioseqexception("Bioseq coordinates in subseq is 1-based, coordinates < 1");
   }
   if (e > (int)length() || e == -1) {
      e = length();
   }

   if (b<=e) {
      bioseq tmp(seq.substr(b-1,e-b+1)); 
      tmp.title = "Subsequence " + itos(b) + "-" + itos(e);
      if (!name.empty()) {
            tmp.title +=  " of " + name;
      }
      tmp.name=name + "sub" + itos(b) + "_" + itos(e);
      return tmp;
   }
   else {
      bioseq tmp(seq.substr(e-1,b-e+1)); 
      tmp.title = "Subsequence " + itos(b) + "-" + itos(e);
      if (!name.empty()) {
         tmp.title += " of " + name ;
      }
      tmp.rvc(); // reverse complement
      tmp.name=name + "sub" + itos(b) + "_" + itos(e) + "rc";
      return tmp;
   }
}

bioseq bioseq::subsequence(int b, int len) const { 
   if (b<0 || (unsigned)b > seq.length()-1) {
      cerr << "seq: " << seq << " b,len: " << b << ',' << len << endl;
      throw bioseqexception("index out of bound in bioseq::subsequence");
   }
   if (len > -1) 
      return bioseq(seq.substr(b, len)); 
   else 
      return bioseq(seq.substr(b)); 
}

// b is 0-based index
bioseq bioseq::subsequenceWithName(int b, int len) const { 
   if (b<0 || (unsigned)b > seq.length()-1) {
      cerr << "seq: " << seq << " b,len: " << b << ',' << len << endl;
      throw bioseqexception("index out of bound in bioseq::subsequence");
   }
   int end=b+len; // for making title
   if (len < 0) end = length();
   string newName = getName() + "sub" + to_string(b+1) + "_" + to_string(end);
   string newTitle = getTitle() + " subsequence " 
         + to_string(b+1) + "-" + to_string(end);
   if (len > -1) 
      return bioseq(newName, seq.substr(b, len), newTitle); 
   else 
      return bioseq(newName, seq.substr(b), newTitle); 
}

int bioseq::locateSubsequence(const string &subseq) const {
   if (subseq.empty()) {
      throw bioseqexception("empty subsequence");
   }
   string::size_type l = seq.find(subseq);
   if (l == string::npos) {
      return -1;
   }
   else {
      return static_cast<int>(l);
   }
}

int bioseq::locateSubsequence(const string &subseq, const size_t b) const {
   if (subseq.empty()) {
      throw bioseqexception("empty subsequence");
   }
   string::size_type l = seq.find(subseq, b);
   if (l == string::npos) {
      return -1;
   }
   else {
      return static_cast<int>(l);
   }
}

int bioseq::locateSubsequence(const string &subseq, const size_t b, const size_t e) const {
   if (subseq.empty()) {
      throw bioseqexception("empty subsequence");
   }
   string::size_type l = seq.find(subseq, b);
   if (l == string::npos) {
      return -1;
   }
   else if (l <= e) {
      return static_cast<int>(l);
   }
   else  return -1;
}

void bioseq::rvc() {
	string tmp;
	tmp.reserve(seq.size()+1);
   string::reverse_iterator it=seq.rbegin();
   char c;
	while (it != seq.rend() ) {
      c=*it;
		if (c == 'a') tmp += 't'; 
      else if (c == 'A') tmp += 'T'; 
		else if (c == 'c') tmp += 'g';
		else if (c == 'C') tmp += 'G';
		else if (c == 'g') tmp += 'c';
		else if (c == 'G') tmp += 'C';
		else if (c == 't') tmp += 'a';
		else if (c == 'T') tmp += 'A';
		else {
			if (isalpha(c)) {
            if (islower(c)) tmp += 'n';
            else tmp += 'N';
         }
			else tmp += c;
		}
		++it;
	}
   title += " reverse complemented";
   seq=tmp;
   if (code != nullptr) { delete[] code; code=nullptr; }
}

map<char,double> bioseq::getFrequency() const {
   map<char, double> frequency;
   for (unsigned int i=0; i<seq.size(); i++) {
      frequency[seq[i]] += 1.0/(double)seq.size();
   }
   return frequency;
}

SequenceType bioseq::guessType() const {
   map<char, double> residueFreq = getFrequency();
   if (residueFreq.size() == 4
         || (residueFreq['A'] + residueFreq['C']  
            + residueFreq['G'] + residueFreq['T'] 
            + residueFreq['N']
            + residueFreq['a'] + residueFreq['c']  
            + residueFreq['g'] + residueFreq['t'] 
            + residueFreq['n']) > 0.95)
      return DNASEQ;
   else if (residueFreq.size() > 12) 
      return PROTEINSEQ;
   else {
      cerr << "Cannot guess sequence type for:\n";
      cerr << *this << endl;
      return UNKNOWN;
   }
}

pair<double,double> bioseq::entropy(const string &seq) {
   map<char, int> frequency;
   map<string, int> freq2aa;
   for (unsigned int i=0; i<seq.size(); i++) {
      frequency[seq[i]]++;
      if (i+1 < seq.length()) 
         freq2aa[seq.substr(i,2)]++;
   }
   map<char, int>::const_iterator it;
   double H=0;
   double f;
   for (it=frequency.begin(); it != frequency.end(); it++) {
      //cout << it->first << " " << it->second << endl;
      f = (double)(it->second)/seq.size();
      H += -f*log(f);
   }
   double H2=0;
   map<string,int>::const_iterator mi;
   //double lensqr = (seq.length()-1)*(seq.length()-1);

   for (mi=freq2aa.begin(); mi != freq2aa.end(); mi++) {
      //f = it->second/lensqr;
      f = (double)mi->second/(seq.length()-1);
      H2 += -f*log(f);
   }
   /* no need to check
   if (H<0 || H2 < 0) {
      cerr << "Entropy cannot be less than zero\n";
      exit(1);
   }
   */
   //return make_pair<double,double>(H,H2);
   return make_pair(H,H2);
}

bioseq& bioseq::operator=(const bioseq &s) {
   if (this != &s) {
      seq=s.seq;
      name=s.name;
      title=s.title;
      if (code != nullptr) {
         delete[] code;
         // don't copy code until it is called with getcode()
         code=nullptr; 
      }
   }
   return *this;
}

bioseq& bioseq::operator=(bioseq &&s) {
   if (this != &s) {
      seq=std::move(s.seq);
      name=std::move(s.name);
      title=std::move(s.title);
      if (code != nullptr) {
         delete[] code;
      }
      code=s.code; 
      s.code=nullptr;
   }
   return *this;
}

bioseq bioseq::operator+(const bioseq &s) const {
   bioseq ss(*this);
   ss += s;
   return ss;
}

bioseq& bioseq::operator=(const string &str) {
   if (str.length()>seq.length())
      seq.resize(str.length());
   seq=str;
   if (code != nullptr) {
      delete[] code;
      code=nullptr;
   }
   return *this;
}

const int* bioseq::getcode() const {
   //cerr << "using bioseq::getcode()\n";
   if (code == nullptr) {
      code = new int[length()+1];
      unsigned int i;
      for (i=0; i<length(); i++) {
         code[i] = aachar2num(seq[i]);
      }
      code[i] = -1; // terminator in case we need for other operation
   }
   return code;
}

void bioseq::toLowerCase() {
   string::iterator i=seq.begin();
   while (i != seq.end()) {
      *i = tolower(*i);
      ++i;
   }
}
void bioseq::toUpperCase() {
   string::iterator i=seq.begin();
   while (i != seq.end()) {
      *i = toupper(*i);
      ++i;
   }
}

void bioseq::show() const {
   cout << name << " ";
   if (!title.empty()) {
      cout << title << " ";
   }
   //printFasta(cout);
   cout << "sequence code\n";
   const int* c = getcode();
   for (unsigned int i=0; i<length(); ++i) {
      cout << "(" << seq[i] << "|" << c[i] << ")";
   }
   cout << endl;
}

void bioseq::appendTitle(const string &aux, const string &sep) { 
   if (title.empty()) {
      title = aux; 
   }
   else {
      title += (sep + aux);
   }
}

void bioseq::assign(string &&s) { 
   seq=std::move(s); 
   if (code != nullptr) { // lazy action
      delete[] code;
      code=nullptr;
   } 
}

void bioseq::setSequence(string &&s) {
   seq=std::move(s); 
   if (code != nullptr) { // lazy action
      delete[] code;
      code=nullptr;
   } 
}

//////////////   DNA Class ////////////////////////////


char DNA::toRevcompBase(const char B) {
   char R;
   switch(B) {
      case 'A': case 'a': R = 'T'; break;
      case 'C': case 'c': R = 'G'; break;
      case 'G': case 'g': R = 'C'; break;
      case 'T': case 't': R = 'A'; break;
      case 'N': case 'n': R = 'N'; break;
      default: R = 'N';
   }
   return R;
}

codon DNA::codontable=codon();

Protein DNA::translate() const {
   string prt;

   for (unsigned int i=0; i<seq.length(); i += 3) {
      prt += codontable[seq.substr(i,3)];
   }
   return Protein(prt);
}

DNA& DNA::operator=(const DNA &s) {
   if (this != &s) {
      bioseq::operator=(s);
   }
   return *this;
}

DNA& DNA::operator=(DNA &&s) {
   if (this != &s) {
      bioseq::operator=(std::move(s));
   }
   return *this;
}

DNA& DNA::operator=(const string &str) {
   bioseq::operator=(str);
   return *this;
}

DNA DNA::subsequence(int b, int len) const { 
   if (b<0 || (unsigned)b > seq.length()-1) {
      cerr << "seq: " << seq << " b,len: " << b << ',' << len << endl;
      throw bioseqexception("index out of bound in bioseq::subsequence");
   }
   if (len > -1) 
      return DNA(seq.substr(b, len)); 
   else 
      return DNA(seq.substr(b)); 
}


// biosequence coordinate
Protein DNA::translate(int begin, int end) const {
   string prt;
   int i=begin-1;
   if (end == -1 || (unsigned)end > seq.length()) 
      end = seq.length();
   while (i< end) {
      prt += codontable[seq.substr(i,3)];
      i += 3;
   }
   return Protein(prt);
}

// code is zeroed on copy
void DNA::revcomp() {
   name += "RC";
   reverseComplementInPlace(seq);
   if (code !=nullptr) {
      delete[] code;
      code=nullptr;
   }
}
// code is zeroed
DNA DNA::revcompCopy() const {
   DNA tmp;
   tmp.name=name + "RC";
   tmp.seq=reverseComplement(seq);
   return tmp;
}

// you are getting a copy not just a pointer
// should return a pointer
const int* DNA::getcode() const {
   //cerr << "Using getcodNuc()\n";
   if (code == nullptr)  {
      unsigned int i;
      code = new int[bioseq::length()+1];
      for (i=0; i<bioseq::length(); i++) {
         code[i]=hashbase(seq[i]);
      }
      code[i] = -1;
   }
   return code;
}

double DNA::GCContent() const {
   int at=0, gc=0;
   for (unsigned int i=0; i<seq.size(); i++) {
      if (seq[i] == 'A' || seq[i] == 'a' || seq[i] == 'T' || seq[i] == 't')
         at++;
      else if (seq[i] == 'G' || seq[i] == 'g' || seq[i] == 'C' || seq[i] == 'c')
         gc++;
   }
   //cerr << "AT: " << at << " GC: " << gc << " GC content: "
    //  << double(gc)/(gc+at) << endl;
   return (double)gc/(gc+at);
}
double DNA::GCContent(long int &A, long int &C, long int &G, long int &T, long int &N) const {
   long int a,c,g,t,n,o;
   a=c=g=t=n=o=0;
   for (unsigned int i=0; i<seq.size(); i++) {
      if (seq[i] == 'A' || seq[i] == 'a') ++a;
      else if (seq[i] == 'T' || seq[i] == 't') ++t;
      else if (seq[i] == 'G' || seq[i] == 'g') ++g;
      else if (seq[i] == 'C' || seq[i] == 'c') ++c;
      else if (seq[i] == 'N' || seq[i] == 'n') ++n;
      else ++o;
   }
   //cerr << "AT: " << at << " GC: " << gc << " GC content: "
    //  << double(gc)/(gc+at) << endl;
   if (o > 0) {
      cerr << o << " bases other than ACGTN\n";
   }
   A += a; C += c; G += g; T += t; N += n;
   return (double)(c+g)/(a+c+g+t);
}

// should return the range in DNA coordinate
// after translation, it could lose information
// if length(DNA) % 3 != 0
pair<unsigned int, unsigned int> longestORFIndex(const Protein &p) {
   string::size_type i, ii, maxb, maxe, len;
   i=0; len=0;
   string ss=p.getSequence();
   ii=ss.find('*', i);
   while (ii != string::npos) {
      if (ii-i > len) {
         len=ii-i;
         maxb=i; maxe=ii+1;
      }
      i=ii+1;
      if (i == ss.length()) break;
      ii=ss.find('*', i);
   }
   if (i < ss.length() && ss.length()-i > len) {
      len=ss.length()-i;
      maxb=i; maxe=ss.length();
   }
   //cout << "max ORF is: from " << maxb << "-" << maxe << endl;
   //cout << ss.substr(maxb, maxe-maxb) << endl;

   return pair<unsigned int, unsigned int>(maxb,maxe);
}

// more efficient version
void DNA::longestORFForward(Protein &p, int &b, int &e) const {
   int maxb, maxe, frame;
   Protein *maxs;
   frame=0;
   Protein p1=translate();
   pair<unsigned int, unsigned int> idx = longestORFIndex(p1);
   maxb=idx.first; maxe=idx.second;
   maxs=&p1;

   Protein p2=translate(2);
   idx=longestORFIndex(p2);
   if (idx.second-idx.first > unsigned(maxe-maxb)) {
      maxb=idx.first; maxe=idx.second;
      frame=1; maxs=&p2;
   }
   Protein p3=translate(3);
   idx=longestORFIndex(p3);
   if (idx.second-idx.first > unsigned(maxe-maxb)) {
      maxb=idx.first; maxe=idx.second;
      frame=2; maxs=&p3;
   }
   p = maxs->subsequence(maxb,maxe-maxb);
   b=frame + 3*maxb;
   if (maxe == (int)maxs->length() && (p[p.length()-1] == '?'
            || p[p.length()-1] == '1' 
            || p[p.length()-1] == '2')) {
      e=length();
   }
   else e=frame + 3*maxe;
}

Protein DNA::longestORFForward(int &b, int &e) const {
   int maxb, maxe, frame;
   Protein *maxs;
   frame=0;
   Protein p1=translate();
   //cout << "frame 0\n" << p1 << endl;
   pair<int,int> idx=longestORFIndex(p1);
   maxb=idx.first; maxe=idx.second;
   maxs=&p1;

   Protein p2=translate(2);
   //cout << "frame 1\n" << p2 << endl;
   idx=longestORFIndex(p2);
   if (idx.second-idx.first > maxe-maxb) {
      maxb=idx.first; maxe=idx.second;
      frame=1; maxs=&p2;
   }
   Protein p3=translate(3);
   //cout << "frame 2\n" << p3 << endl;
   idx=longestORFIndex(p3);
   if (idx.second-idx.first > maxe-maxb) {
      maxb=idx.first; maxe=idx.second;
      frame=2; maxs=&p3;
   }
   Protein lp = maxs->subsequence(maxb,maxe-maxb);
   //cout << "longest of three frames: " << frame << endl << lp << endl;
   b=frame + 3*maxb;
   if (lp[lp.length()-1] == '?') {
      e=length();
   }
   else e=frame + 3*maxe;
   return lp;
}

void longestORFPlus(const string &rna, string &pep, int &b, int &e) {
   int f;
   longestORFPlus(rna,pep,b,e,f);
   if (f != b % 3) {
      cerr << rna << endl;
      cerr << "begin: " << b << " b%3 " << b%3 << endl;
      cerr << "frame: " << f << " should be derivable from begin!\n";
      exit(1);
   }
}

void longestORFPlus(const string &rna, string &pep, int &b, int &e, int &f) {
   vector<Interval> candidates;
   //static const int cutoff=5; // AA length cutoff
   int i;

   string peps[3];
   string::size_type m,s; // start and stop
   for (i=0; i<3; i++) {
      translate(peps[i], rna, i+1, -1); // translate use 1- index
      m=peps[i].find('M');
      if (m == string::npos) { // no M in the whole seq
         if ((s=peps[i].find('*')) == string::npos)
            m=0;
         else { // no start but has stop
            candidates.push_back(Interval(i, P2Rindex(s,i)+2));
            continue;
         }
      }
      else {// found m
         if (m>30 && peps[i].substr(0,m).find('*') == string::npos) 
         { m=0; }
      }
      while (m != string::npos ) {
         s=peps[i].find('*', m+1);
         if (s != string::npos) {
            candidates.push_back(P2Rindex(Interval(m,s), i));
         }
         else {
            candidates.push_back(Interval(P2Rindex(m,i),rna.length()-1));
            break;
         }
         m=peps[i].find('M', s+1);
      }
   }
   if (candidates.empty()) {
      cerr << "no ORF found for RNA:\n" << rna << endl;
      exit(1);
   }
   Interval maxo;
   for (i=0; (unsigned)i < candidates.size(); ++i) {
      if (candidates[i].length() > maxo.length()) {
         maxo=candidates[i];
      }
   }
   // max protein is pB=n/3, at frame n%3 start
   // end if n=rna.length()-1, pE=pep.end, else pB=n/3
   if ((unsigned)maxo.end() == rna.length()-1) {
      pep=peps[maxo.begin()%3].substr(maxo.begin()/3); 
   }
   else 
      pep=peps[maxo.begin()%3].substr(maxo.begin()/3, maxo.length()/3); 
   b=maxo.begin();
   e=maxo.end();
   f=b%3;
}

void longestORFPlusSuffix(const string &rna, string &pep, int &b, int &e) {
   unsigned int suff, midlen=0; 
   int suffframe, midframe;
   pair<int,int> mid;
   suffframe=midframe=-1;

   // pre is the stop codon position, 
   // mid is the middle orf, suff is the start 
   // codon of the tail

   string peps[3];
   string::size_type m,s; // start and stop
   suff=rna.length()/3;
   for (int i=0; i<3; i++) {
      translate(peps[i], rna, i+1, -1);
      m = peps[i].find('M');
      while (m != string::npos) {
         s=peps[i].find('*', m+1);
         if (s == string::npos) { // tail ORF
            if (m < suff) {
               suff = m;
               suffframe=i;
            }
            break;
         }
         else {
            if (s-m > midlen) {
               midlen=s-m;
               mid=make_pair(m,s);
               midframe=i;
            }
         }
         ++s;
         if (s == peps[i].length()) break;
         m = peps[i].find('M', s);
      }
   }
   unsigned int sufflen=0;
   if (suffframe != -1) sufflen=peps[suffframe].length()-suff;
   if (midlen > sufflen) {
      b=midframe + 3*mid.first;
      e=midframe + 3*mid.second+2;
      pep=peps[midframe].substr(mid.first, mid.second-mid.first+1);
   }
   else if (sufflen > midlen) {
      b=suffframe + 3*suff;
      e=rna.length()-1;
      pep=peps[suffframe].substr(suff);
   }
}

void longestORFPlusPrefix(const string &rna, string &pep, int &b, int &e) {
   unsigned int pre=0, midlen=0; 
   int preframe, midframe;
   pair<int,int> mid;
   preframe=midframe=-1;

   // pre is the stop codon position, 
   // mid is the middle orf, suff is the start 
   // codon of the tail

   string peps[3];
   string::size_type m,s; // start and stop
   for (int i=0; i<3; i++) {
      translate(peps[i], rna, i+1, -1);
      s=peps[i].find('*');
      if (s == string::npos) {
         cerr << "whole pepseq no stop, never seen before\n";
         exit(1);
      }
      if (s > pre)  {
         pre=s;
         preframe=i;
      }
      m = peps[i].find('M', s+1);
      while (m != string::npos) {
         s=peps[i].find('*', m+1);
         if (s == string::npos) { // tail ORF
            break;
         }
         else {
            if (s-m > midlen) {
               midlen=s-m;
               mid=make_pair(m,s);
               midframe=i;
            }
         }
         ++s;
         if (s == peps[i].length()) break;
         m = peps[i].find('M', s);
      }
   }
   if (midlen > pre) {
      b=midframe + 3*mid.first;
      e=midframe + 3*mid.second+2;
      pep=peps[midframe].substr(mid.first, mid.second-mid.first+1);
   }
   else if (pre > midlen) {
      b=0;
      e=preframe + 3*pre+2;
      pep=peps[preframe].substr(0,pre+1);
   }
}

Protein DNA::longestORFReverse(int &b, int &e) const {
   DNA tmp=*this;
   tmp.revcomp();
   Protein p=tmp.longestORFForward(b,e);
   //cout << "length: " << length() << endl;
   b=length()-b;
   e=length()-e;
   return p;
}

void DNA::longestORFReverse(Protein &p, int &b, int &e) const {
   DNA tmp=*this;
   tmp.revcomp();
   tmp.longestORFForward(p, b,e);
   //cout << "length: " << length() << endl;
   b=length()-b;
   e=length()-e;
}

Protein DNA::longestORF(char strand, int &b, int &e) const {
   if (strand == '+') return longestORFForward(b,e);
   else return longestORFReverse(b,e);
}

Protein DNA::longestORF(int &b, int &e) const {
   Protein p1=longestORFForward(b,e);
   int bb, ee;
   Protein p2=longestORFReverse(bb,ee);
   if (bb-ee <= e-b) {
      //cout << b << "  " << e << " longer\n";
      return p1;
   }
   else {
      //cout << bb << "  " << ee << " longer\n";
      b=bb; e=ee;
      return p2;
   }
}

void DNA::longestORF(Protein &p, int &b, int &e) const {
   longestORFForward(p, b,e);
   int bb, ee;
   Protein p2;
   longestORFReverse(p2, bb, ee);
   if (bb-ee <= e-b) {
      cout << b << "  " << e << " longer\n";
   }
   else {
      cout << bb << "  " << ee << " longer\n";
      b=bb; e=ee;
      p=p2;
   }
}
bool DNA::noAmbiguous() const {
   for (size_t i=0; i<seq.size(); ++i) {
      if (toupper(seq[i]) != 'A'
            && toupper(seq[i]) != 'C'
            && toupper(seq[i]) != 'G'
            && toupper(seq[i]) != 'T')
         return false;
   }
   return true;
}
bool DNA::ambiguous() const {
   for (size_t i=0; i<seq.size(); ++i) {
      if (toupper(seq[i]) != 'A'
            && toupper(seq[i]) != 'C'
            && toupper(seq[i]) != 'G'
            && toupper(seq[i]) != 'T')
         return true;
   }
   return false;
} 

bool DNA::read(istream &ins) {
   if (ins.eof()) return false;
   string::size_type i;
   //clear(); // make sure the sequence starts with empty
   getline(ins,title);
   if (title[0] != '>' || title.length() < 2) {
      cerr << "not in proper fasta format, first line of DNA file: "
         << title << endl;
      throw bioseqexception("fasta format problem");
   }
   //title=title.substr(1); // get rid of >
   if ((i=title.find(' ')) != string::npos) {
      name=title.substr(1,i-1);
      title=title.substr(i+1);
      //if (title[0] == ' ') { title=title.substr(1); }
   }
   else {
      name=title;
      title.clear();
   }
   if (code != nullptr) {
      delete[] code;
      code=nullptr;
   }
   seq.clear();
   string line; 
   getline(ins,line);
   while (!ins.eof()) {
      if (line.empty()) {
         getline(ins,line);
         continue;
      }
      while (!isprint(line.back())) { // remove invisible variable
         line.resize(line.length()-1);
      }
#ifdef VERIFY_RESIDUE
      if (line.find_first_not_of("ACGTURYSWKMBDHVN-acgturyswkmbdhvn") != string::npos) 
      {
         throw bioseqexception("None nucleotide symbol in line: " + line + "|");
      }
#endif
      if (!line.empty() && 
            (line.find('U') != string::npos || line.find('u') != string::npos)) {
         for (size_t j=0; j < line.size(); ++j) {
            if (line[j] == 'U') line[j] = 'T';
            else if (line[j] == 'u') line[j] = 't';
         }
      }
      seq += line;
      if (ins.peek() == '>') return true;
      getline(ins,line);
   }
   if (seq.empty()) {
      cerr << "Empty sequence something is wrong!\n";
      throw bioseqexception("empty sequence from sequence file");
   }
   if (!ins.eof())  {
      ins.peek();
   }
   return true;
}


bool longestNoStartORFPlus(const string &rna, 
      pair<int,int> &nterm, string &npep, pair<int,int> &full, string &pep)
{
   // ns: no start
   // E end, B begin, L length, F frame
   int i, Bfull, Efull, Ffull, Fns;
   unsigned int Ens;
   Ens=Bfull=Efull=0;
   Ffull=Fns=-1;
   int NSNS_frame=-1;

   string peps[3];
   string::size_type m,s; // start and stop
   for (i=0; i<3; i++) {
      translate(peps[i], rna, i+1, -1);
      s=peps[i].find('*'); m=peps[i].find('M');
      if (s == string::npos) { // usually short sequences!
         NSNS_frame=i; // with M or not
         continue;
      }
      else { // has stop
         if (m == string::npos) { // no start, NoStart ORF
            // --------*----- 
            if (s > Ens)  { Ens=s; Fns=i; }
            continue;
         }
         else { // has start
            if (m>s && s > Ens)  { Ens=s; Fns=i; }
         }
      }
      while (m != string::npos) {
         s=peps[i].find('*', m+1);
         if (s == string::npos) { // --M--- Not interested
            break;
         }
         else {
            if (s-m > (unsigned)(Efull - Bfull)) {
               Bfull=m; Efull=s; Ffull=i;
            }
         }
         ++s;
         if (s >= peps[i].length()-1) break;
         m = peps[i].find('M', s);
      }
   }
   if (Ens == 0 && Efull == 0) {
      if (NSNS_frame != -1) {
         pep=peps[NSNS_frame];
         full=make_pair(NSNS_frame, rna.length()-1);
      }
      return false;
   }

   if (Efull > Bfull) {
      full=make_pair(Ffull + 3*Bfull, Ffull + 3*Efull+2);
      pep=peps[Ffull].substr(Bfull, Efull - Bfull + 1);
   }
   else {
      full=make_pair(0,0);
      pep="";
   }
   if (Ens > 0) {
      nterm=make_pair(Fns, Fns+3*Ens+2);
      npep=peps[Fns].substr(0,Ens+1);
   }
   else {
      nterm=make_pair(0,0);
      npep="";
   }
   return true;
}

// pick full and no stop ORF
bool longestNoStopORFPlus(const string &rna, 
      pair<int,int> &cterm, string &cpep, pair<int,int> &full, string &pep)
{
   // ns: no stop
   // E end, B begin, L length, F frame
   int i, Bfull, Efull, Ffull, Fns, NSNS_frame;
   Bfull=Efull=0;
   Ffull=Fns=NSNS_frame=-1;
   unsigned int Bns=rna.length()/3;

   string peps[3];
   string::size_type m,s; // start and stop
   for (i=0; i<3; i++) {
      translate(peps[i], rna, i+1, -1);
      s=peps[i].find('*'); m=peps[i].find('M');
      if (s == string::npos) { // usually short sequences!
         if (m == string::npos) {
            NSNS_frame=i; // ------- No (M & *)
         }
         else { // ---M-----
            if (m < Bns) {
               Bns = m; Fns=i;
            }
         }
         continue;
      }
      else { // has stop
         if (m == string::npos) { // no start, NoStart ORF
            // --------*----- 
            continue;
         }
         //else { // has start ---*---M--?-. do nothing
         //leave to the while loop to process
      }
      while (m != string::npos) {
         s=peps[i].find('*', m+1);
         if (s == string::npos) { // --M--- collect end case
            if (m < Bns) {
               Bns=m; Fns=i;
            }
            break;
         }
         else {
            if (s-m > unsigned(Efull - Bfull)) {
               Bfull=m; Efull=s; Ffull=i;
            }
         }
         ++s;
         if (s >= peps[i].length()-1) break;
         m = peps[i].find('M', s);
      }
   }
   if (Bns == rna.length()/3 && Efull == 0) {
      if (NSNS_frame != -1) {
         pep=peps[NSNS_frame];
         full=make_pair(NSNS_frame, rna.length()-1);
      }
      return false;
   }

   if (Efull > Bfull) {
      full=make_pair(Ffull + 3*Bfull, Ffull + 3*Efull+2);
      pep=peps[Ffull].substr(Bfull, Efull - Bfull + 1);
   }
   else {
      full=make_pair(0,0);
      pep="";
   }
   if (Bns < rna.length()/3) {
      cterm=make_pair(Fns + 3*Bns, rna.length()-1);
      cpep=peps[Fns].substr(Bns);
   }
   else {
      cterm=make_pair(0,0);
      cpep="";
   }
   return true;
}

///////////// Protein Class /////////////////////////
//
/* Amino acid symbols
One-letter Three-letter Amino acid
----------------------------------------------------
A  Ala   alanine
B  Asx   aspartic acid or asparagine
C  Cys   cysteine
D  Asp   aspartic acid
E  Glu   glutamic acid
F  Phe   phenylalanine
G  Gly   glycine
H  His   histidine
I  Ile   isoleucine
J*  NUL   not defined***
K  Lys   lysine
L  Leu   leucine
M  Met   methionine
N  Asn   asparagine
O*  NUL   not defined ***
P  Pro   proline
Q  Gln   glutamine
R  Arg   arginine
S  Ser   serine
T  Thr   threonine
U* Sec   selenocysteine
V  Val   valine
W  Trp   tryptophan
X**   Xaa   unknown or 'other' amino acid
Y  Tyr   tyrosine
Z  Glx   glutamic acid or glutamine (or substances such as
4-carboxyglutamic acid and 5-oxoproline that
yield glutamic acid on acid hydrolysis of peptides)
*/

const char * Protein::symbols[]={
"Ala", "alanine", 
"Asx", "aspartic acid or asparagine",
"Cys", "cysteine", 
"Asp", "aspartic acid",
"Glu", "glutamic acid", 
"Phe", "phenylalanine",
"Gly", "glycine", 
"His", "histidine", 
"Ile", "isoleucine",
"NUL", "not defined",
"Lys", "lysine", 
"Leu", "leucine", 
"Met", "methionine",
"Asn", "asparagine", 
"NUL", "not defined",
"Pro", "proline", 
"Gln", "glutamine",
"Arg", "arginine", 
"Ser", "serine", 
"Thr", "threonine",
"Sec", "selenocysteine", 
"Val", "valine", 
"Trp", "tryptophan",
"Xaa", "unknown or 'other' amino acid", 
"Tyr", "tyrosine",
"Glx", "glutamic acid or glutamine (or substances such as 4-carboxyglutamic acid and 5-oxoproline that yield glutamic acid on acid hydrolysis of peptides)",
"Stp", "stop codon"};

const char* Protein::oneLetterSymbols="ABCDEFGHIKLMNPQRSTUVWXYZ*";

/* 
 * Release 49.0 of 07-Feb-2006 of UniProtKB/Swiss-Prot contains 207'132
 * sequence entries,
 * comprising 75'438'310 amino acids abstracted from 139'151 references.
 *
 * 12'815 sequences have been added since release 48, the sequence data of 
 * 991 existing 
 * entries has been updated and the annotations of all entries have been
 * revised. 
 * This represents an increase of 7%.
 * percent
Ala (A) 7.83   Gln (Q) 3.95   Leu (L) 9.64   Ser (S) 6.86
Arg (R) 5.35   Glu (E) 6.64   Lys (K) 5.93   Thr (T) 5.42
Asn (N) 4.18   Gly (G) 6.93   Met (M) 2.38   Trp (W) 1.15
Asp (D) 5.32   His (H) 2.29   Phe (F) 4.00   Tyr (Y) 3.06
Cys (C) 1.52   Ile (I) 5.91   Pro (P) 4.83   Val (V) 6.71

map<char, double> AAFrequency = map<char,double>;
char AAlist[]={'A', 'Q', 'L', 'S', 'R', 'E', 'K', 'T',
         'N', 'G', 'M', 'W', 'D', 'H', 'F', 'Y',
         'C', 'I', 'P', 'V'}; 
double AAfreq[]= {0.0783, 0.0395, 0.0964, 0.0686, 0.0535,
        0.0664, 0.0593, 0.0542, 0.0418, 0.0693, 0.0238,
        0.0115, 0.0532, 0.0229, 0.0400, 0.0306, 0.0152,
        0.0591, 0.0483, 0.0671};

*/

const double Protein::aafreq[] = { 
   0.0783, 0.0000, 0.0152, 0.0532, 0.0664,
   0.0400, 0.0693, 0.0229, 0.0591, 0.0000,
   0.0593, 0.0964, 0.0238, 0.0418, 0.0000,
   0.0483, 0.0395, 0.0535, 0.0686, 0.0542, 
   0.0000, 0.0671, 0.0115, 0.0000, 0.0306, 0.0000};

Protein& Protein::operator=(const Protein &p) {
   if (this != &p) {
      bioseq::operator=(p);
      for (int i=0; i<numcodes; i++)
         codefreq[i]=p.codefreq[i];
   }
   return *this;
}

Protein& Protein::operator=(Protein &&p) {
   if (this != &p) {
      bioseq::operator=(std::move(p));
      codefreq=std::move(p.codefreq);
      //for (size_t i=0; i<numcodes; i++)
      //   codefreq[i]=p.codefreq[i];
   }
   return *this;
}

Protein& Protein::operator=(const string &str) {
   bioseq::operator=(str);
   codefreq[0]=-1;
   return *this;
}

double Protein::relativeEntropy() const {
   double tmp=0;
   map<char, double> freq=getFrequency();
   double f[26]={0.0};
   map<char, double>::const_iterator it;
   for (it=freq.begin(); it != freq.end(); it++) {
      //cout << it->first << " --> " << it->second << endl;
      f[toupper(it->first)-'A']=it->second;
   }

   for (int i=0; i<26; i++) {
      if (aafreq[i] <= 0 || f[i] <= 0) continue;
      //cout << char('A'+i) << " " << f[i] << " " 
       //  << aafreq[i] << endl;
      tmp += f[i]*log(f[i]/aafreq[i]);
   }
   return tmp;
}

// more efficient version
//const double* Protein::getCodeFrequence() {
array<double, 27>& Protein::getCodeFrequence() {
   const int *c = getcode();
   //while (*c > -1) cerr << *c << " ";
   //cerr << endl;
   if (codefreq[0] < 0.0) {
      for(int i=0; i<27; i++) codefreq[i]=0.0;
      while (*c != -1) {
         codefreq[*c] += 1.00/(double)length();
         ++c;
      }
   }
   return codefreq;
}

double Protein::relativeEntropyUniform() const {
   double tmp=0;
   map<char, double> freq=getFrequency();
   double f[26]={0.0};
   map<char, double>::const_iterator it;
   for (it=freq.begin(); it != freq.end(); it++) {
      f[toupper(it->first)-'A']=it->second;
   }
   codon cc;
   map<char,double> codonf=cc.getAAUniformFrequency();
   double uf[26]={0.0};
   for (it=codonf.begin(); it != codonf.end(); it++) {
      uf[it->first-'A']=it->second;
   }

   for (int i=0; i<26; i++) {
      if (uf[i] <= 0 || f[i] <= 0) continue;
      tmp += f[i]*log(f[i]/uf[i]);
   }
   return tmp;
}

bool Protein::hasInternalStop() const {
   string::size_type i = seq.find("*");
   if (i == string::npos || i != seq.length()-1) {
      return false;
   }
   else return true;
}

int Protein::countInternalStops() {
   string::size_type i=0, j;
   int count=0;
   while (i < seq.size()-1 && (j=seq.find("*", i)) != string::npos) {
      ++count;
      i=j+1;
   }
   return count;
}

bool loadFastaIntoMap(const string &file, map<string,string> &store) {
   ifstream IN(file.c_str());
   if (IN.fail()) {
      cerr << "failed to open fasta file " << file << endl;
      return false;
      //exit(1);
   }
   string line, id, seq;
   getline(IN, line);
   string::size_type i;
   while (!IN.eof() && line[0] == '>') {
      i=line.find(' ');
      if (i != string::npos) {
         id=line.substr(1,i-1);
      }
      else id=line.substr(1);
      getline(IN,line);
      if (!seq.empty()) seq.clear();
      while (!IN.eof() && line[0] != '>') {
         seq += line;
         getline(IN,line);
      }
      store[id]=seq;
   }
   return true;
}


///// mRNA class //////////////

mRNA::mRNA(const string &s) : DNA(s), prt() {
   longestORF(prt, cdsb, cdse);
}

//////////// DNAQual ////////////

int DNAQual::defaultQ=60;

DNAQual::DNAQual(const DNAQual &s) : DNA(s), qual(0), rc(0) {
   if (s.qual != 0) {
      qual = new int[length()];
      for (size_t i=0; i<length(); ++i) {
         qual[i]=s.qual[i];
      }
   }
}

DNAQual::DNAQual(const string &n, const string &s, int *qs)
   : DNA(n,s), qual(new int[length()]), rc(0) {
      for (size_t i=0; i<length(); ++i) {
         qual[i]=qs[i];
      }
}

DNAQual::DNAQual(const string &n, const string &s, const unsigned char *qs)
   : DNA(n,s), qual(new int[length()]), rc(0) {
      for (size_t i=0; i<length(); ++i) {
         qual[i]=int(qs[i]);
      }
}

DNAQual::DNAQual(const string &n, const string &s, const vector<int> &qs)
   : DNA(n,s), qual(new int[length()]), rc(0) {
      for (size_t i=0; i<length(); ++i) {
         qual[i]=qs[i];
      }
}

DNAQual& DNAQual::operator=(const DNAQual &dq) {
   if (this != &dq) {
      // first take care of quality
      if (dq.qual == 0) { // no qualty from source
         if (qual != 0) {
            delete[] qual;
            qual=0;
         }
      }
      else { // need to have enough space for quality
         if  (qual != 0) {
            if (length() < dq.length()) {
               delete[] qual;
               qual = new int[dq.length()];
            }
         }
         else qual = new int[dq.length()];
      }
      DNA::operator=(dq);
      if (qual != 0) {
         for (size_t i=0; i<dq.length(); ++i) {
            qual[i] = dq.qual[i];
         }
      }
      // never copy rc object, keep rc empty for new objects
      if (rc != 0) { 
         delete rc;
         rc = 0;
      }
   }
   return *this;
}

DNAQual& DNAQual::operator=(DNAQual &&dq) {
   if (this != &dq) {
      DNA::operator=(std::move(dq));
      if (qual != 0) delete[] qual;
      qual=dq.qual;
      dq.qual=0;
      rc=dq.rc;
      dq.rc=0;
   }
   return *this;
}

DNAQual& DNAQual::operator=(const string &str) { 
   setSequence(str);
   if (qual != 0) {
      delete[] qual;
      qual = 0;
   }
   if (rc != 0) {
      delete rc;
      rc=0;
   }
   return *this;
}

DNAQual& DNAQual::operator=(const DNA &dna) { 
   if (this != &dna) {
      DNA::operator=(dna); 
      if (qual != 0) {
         delete[] qual;
         qual = 0;
      }
      if (rc != 0) {
         delete rc;
         rc = 0;
      }
   }
   return *this;
}

DNAQual::~DNAQual() { 
   if (qual != 0) {
      delete[] qual; 
      qual=0; 
   }
   if (rc != 0) {
      delete rc; 
   }
}

DNAQual DNAQual::subseq(int b, int e) const {
   DNAQual tmp(DNA::subseq(b,e));
   int i=b-1;
   int x=0;
   tmp.qual=new int[e-b+1]; 
   while (i<e) {
      tmp.qual[x]=qual[i];
      ++x; ++i;
   }
   return tmp;
}

DNAQual DNAQual::subsequenceWithName(int b, int len) const {
   int endPos = len<0? length() : b+len;
   DNAQual tmp = this->subseq(b+1, endPos);
   tmp.setName(tmp.name + "_" + itos(b) + "_" + itos(b+len));
   if (tmp.hasTitle()) {
      tmp.setTitle(tmp.getTitle() + " subsequence from " + itos(b) 
            + " to " + itos(b+len));
   }
   else {
      tmp.setTitle(" subsequence from " + itos(b) + " to " + itos(b+len));
   }
   //if (rc != 0) {
   //   delete rc;
   //   rc=0;
   //}
   return tmp;
}

void DNAQual::reverseQuality() {
   if (qual == 0) return;
   size_t i, L;
   L = length() - 1;
   for (i=0; i<length()/2; ++i) 
      swap(qual[i], qual[L-i]);
}
// the code is zerod on copy construction and copy assignment
// so it has to be dealt with separately
void DNAQual::revcomp() {
   if (rc == 0) {
      rc=new DNAQual(*this); // rc is a copy of the current object
      DNA::revcomp(); // rc current object sequence
      reverseQuality();
   }
   else {
      /* slower working
      *this=*rc;
      delete rc;
      rc=0;
      */
      //cout << "Inside revcomp(): swaping rc and this\n";
      //cout << "before swap\n" << *this << endl;
      int *rcCode=rc->code;
      rc->code=nullptr;
      int *thisCode=code;
      code=nullptr;
      DNAQual* rcptr=rc;
      this->rc=nullptr;
      DNAQual* thisPtr=new DNAQual(*this); // can this be moved?
      operator=(std::move(*rcptr)); // rc is this object now
      rc = thisPtr;
      //rc->rc=0; // break circular reference.
      code=rcCode;
      rc->code=thisCode;
      //cout << "after swap\n" << *this << endl;
   }

}

DNAQual DNAQual::revcompCopy() const {
   if (rc != 0) {
      return DNAQual(*rc);
   } // rc has value
   DNAQual tmp(*this); // tmp.rc=0
   tmp.revcomp();
   return tmp;
}

// more efficient version
DNAQual* DNAQual::getRevcomp() {
   if (rc == 0) {
      rc=new DNAQual(*this); // rc is a copy of the current object
      rc->DNA::revcomp(); // rc current object sequence
      reverseQuality();
   }
   return rc;
}

void DNAQual::setQuality(int *qs) {
   for (size_t i=0; i<length(); ++i) {
      qual[i]=qs[i];
   }
}

// will not print the RC information
// this method will be used for input/output to file
ostream& DNAQual::print(ostream &ous, int width) const {
   bioseq::print(ous,width);
   if (qual != 0) {
      for (size_t i=0; i<length(); ++i) {
         if ((i+1)%width == 0) ous << qual[i] << endl;
         else {
            if (i == length()-1) 
               ous << qual[i] << endl;
            else
               ous << qual[i] << ' ';
         }
      }
      //if (i != length()) ous << endl;
   }
   else { // use default score
      for (size_t i=0; i<length(); ++i) {
         if ((i+1)%width == 0) ous << defaultQ << endl;
         else {
            if (i == length()-1) 
               ous << defaultQ << endl;
            else
               ous << defaultQ << ' ';
         }
      }
   }
   /*
   if (rc == 0) {
      ous << "no RC object cached yet\n";
   }
   else {
      ous << "has a cached rc object\n"; 
      ous << *rc << endl;
   }
   */
   return ous;
}

/////// DNAQualCount ////////////
//
DNAQualCount::DNAQualCount(DNA &&dna, int n, const vector<int> &qs)
   : DNAQual(std::move(dna)), cnt(n) 
{
   qual = new int[length()];
   for (size_t i=0; i<length(); ++i) {
      qual[i]=qs[i];
   }
}

DNAQualCount::DNAQualCount(const DNA &dna, int n, const vector<int> &qs)
   : DNAQual(dna), cnt(n) 
{
   qual = new int[length()];
   for (size_t i=0; i<length(); ++i) {
      qual[i]=qs[i];
   }
}


DNAQualCount& DNAQualCount::operator=(const DNAQualCount &dq) {
   if (this != &dq) {
      DNAQual::operator=(dq);
      cnt =dq.cnt;
   }
   return *this;
}

DNAQualCount& DNAQualCount::operator=(const DNAQual &dq) {
   if (this != &dq) {
      DNAQual::operator=(dq);
      cnt = 1;
   }
   return *this;
}

DNAQualCount& DNAQualCount::operator=(DNAQualCount &&dq) {
   if (this != &dq) {
      DNAQual::operator=(std::move(dq));
      cnt=dq.cnt;
   }
   return *this;
}

DNAQualCount DNAQualCount::subseq(int b, int e) const {
   DNAQualCount tmp(DNAQual::subseq(b,e));
   tmp.cnt = cnt;
   return tmp;
}

DNAQualCount DNAQualCount::subsequenceWithName(int b, int len) const {
   DNAQualCount tmp(DNAQual::subsequenceWithName(b, len));
   tmp.cnt = cnt;
   return tmp;
}

DNAQualCount DNAQualCount::revcompCopy() const {
   DNAQualCount tmp(*this);
   // this is base class (DNAQual) method
   tmp.revcomp();
   return tmp;
}

DNAQualCount* DNAQualCount::getRevcomp() {
   if (rc == 0) {
      rc = new DNAQualCount(*this);
      rc->DNA::revcomp();
      rc->reverseQuality();
   }
   return (DNAQualCount*)rc;
}

ostream& DNAQualCount::print(ostream &ous, int width) const {
   DNAQual::print(ous,width); // sequence, then quality number
   ous << "cnt=" << cnt << "\n";
   return ous;
}

istream& operator>>(istream& ins, DNAQualCount &sq) {
   //if (ins.eof()) {
   //   cerr << "end of file\n";
   //   sq.clear();
   //   return ins;
   //}
   string line;
   getline(ins, line);
   //if (ins.eof()) {
   //   cerr << "end of file after attempt to read\n";
   //   sq.clear();
   //   return ins;
   //}
   //cout << line << endl;
   if (line[0] != '>') {
      throw runtime_error("DNAQualCount object text representaton wrong");
   }
   size_t i=line.find(' ');
   if (i == string::npos) {
      sq.setName(line.substr(1));
   }
   else {
      sq.setName(line.substr(1, i-1));
      sq.setTitle(line.substr(i+1));
      if (sq.title[0] == ' ') {
         sq.title = sq.title.substr(1);
      }
   }
   string tmpseq;
   // read sequence part
   getline(ins, line);
   while (!ins.eof() && !isdigit(line[0]) 
         && line.find_first_of("ACGTacgt") != string::npos) {
      tmpseq += line;
      getline(ins, line);
   }
   //cout << tmpseq << endl;
   if (sq.length() < tmpseq.length()) {
      if (sq.qual != 0) {
         delete[] sq.qual;
      }
      sq.qual=new int[tmpseq.length()];
   }
   sq.setSequence(tmpseq);
   //cout << line << endl; // first number line
   i=0;
   while (!ins.eof() && line.substr(0,4) != "cnt=") {
      if (line[line.size()-1] == ' ') {
         line = line.substr(0,line.size()-1);
      }
      vector<string> qv=split(line, ' ');
      //cout << line << endl;
      //cout << qv.size() << endl;
      for (size_t x=0; x<qv.size(); ++x) {
         sq.qual[i++]=stoi(qv[x]);
      }
      getline(ins, line);
   }
   //cout << "should be the count line: " << line << endl;
   sq.setCount(stoi(line.substr(4)));
   ins.peek();
   return ins;
}

ostream& DNAQualCount::printFastaWithHeader(ostream &ous, unsigned int width) const  
{
   string titleTxt= to_string(getCount());
   if (getCount() == 1) titleTxt += " copy";
   else titleTxt += " copies";
   if (!title.empty()) { 
      titleTxt += (", " + title); 
   }
   ous << ">" << getName() << " " << titleTxt << endl;
   printFasta(ous,width);
   return ous;
}
}
