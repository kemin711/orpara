#include "braboualn.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <stddev.h>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <map>
#include <bioseq.h>

using namespace std;
using namespace orpara;

/**
 * For parameters of the program.
 */
struct Progparam {
   public:
      Progparam() : m(5), s(-3), go(-23), ge(-2), ic(0.72) { }
      Progparam(int match, int mismatch, int gapo, int gape, double idcut) 
         : m(match), s(mismatch), go(gapo), ge(gape), ic(idcut) { validate(); }
      void validate() {
         if (m <= 0) {
            cerr << "match score must be greater than zero\n";
            exit(1);
         }
         if (s > 0) {
            cerr << "mis match score should be a negative number\n";
            s = -s;
         }
         if (go > 0) {
            cerr << "gap open score should be a negative number\n";
            go = -go;
         }
         if (ge > 0) {
            cerr << "gap extension score should be a negative number\n";
            ge = -ge;
         }
      }
      int m, s, go, ge;
      double ic;
};

void validateScores(int &mis, int &go, int &ge);
pair<string,string> readSequence(const string &fname);
void usage(const Progparam &param) {
   cerr << "bbaln -1 seq1.fasta -2 seq2.fasta -o outputfile\n"
      << " or bbaln -o outfile seq1.fasta seq2.fasta\n"
      << "Options:\n"
      << "   -m or --match match score. Default " << param.m << "\n"
      << "   -s or --mismatch mismatch score. Default " << param.s << "\n"
      << "   -g or --gap-open gap open score. Default " << param.go << "\n"
      << "   -e or --gap-extend gap extension score. Default " << param.ge << "\n"
      << "   -o or --outfile output file name\n"
      << "   -i or --identity-cut identity cutoff. Default " << param.ic << "\n"
      << "   -h, ?, or --help to show help message\n";
   exit(1);
}

string makeOutFileName(const string &file1, const string &file2);
string extractSeqid(const string &txt);
vector<pair<string,string> > readMultipleSequences(const string &fname);
int selfAlign(const string &infile, const Progparam &param);
vector<DNA> slurpDNA(const string &fname);

/**
 * Actual alignment program for user to use the algorithm.
 */
int main(int argc, char* argv[]) {
   string seq1file, seq2file, outfile;
   //int match=4, mismatch=-3, gapo=-5, gape=-2;
   Progparam param(9, -7, -37, -1, 0.61);
   bool SELF_ALIGN=false;
   int i=1;
   while (i < argc) {
      if (!strcmp(argv[i], "-1")) seq1file=argv[++i];
      else if (!strcmp(argv[i], "-2")) seq2file=argv[++i];
      else if (!strcmp(argv[i], "--mismatch") || !strcmp(argv[i], "-s")) 
         param.s=atoi(argv[++i]);
      else if (!strcmp(argv[i], "--match") || !strcmp(argv[i], "-m")) 
         param.m=atoi(argv[++i]);
      else if (!strcmp(argv[i], "--gap-open") || !strcmp(argv[i], "-g")) 
         param.go=atoi(argv[++i]);
      else if (!strcmp(argv[i], "--gap-extend") || !strcmp(argv[i], "-e")) 
         param.ge=atoi(argv[++i]);
      else if (!strcmp(argv[i], "--outfile") || !strcmp(argv[i], "-o")) 
         outfile=string(argv[++i]);
      else if (!strcmp(argv[i], "--identity-cut") || !strcmp(argv[i], "-i")) 
         param.ic=atof(argv[++i]);
      else if (!strcmp(argv[i], "--self")) {
         SELF_ALIGN=true;
         seq1file=string(argv[++i]);
      }
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")
            || !strcmp(argv[i], "?")) 
         usage(param);
      else {
         if (argv[i][0] == '-') {
            cerr << "wrong optins\n";
            usage(param);
         }
         if (i+1 < argc) {
            seq1file=string(argv[i++]);
            seq2file=string(argv[i]);
         }
         else {
            usage(param);
         }
      }
      ++i;
   }
   if (SELF_ALIGN) { // align all sequences inside the file
      return selfAlign(seq1file, param);
   }
   if (seq1file.empty() || seq2file.empty()) { usage(param); }
   if (outfile.empty()) {
      outfile = makeOutFileName(seq1file, seq2file);
   }
   pair<string, string> seq1 = readSequence(seq1file);
   pair<string, string> seq2 = readSequence(seq2file);

   time_t timer1, timer2;
   time(&timer1);
   Braboualn aligner(seq1.second, seq2.second, param.m, param.s, param.go, param.ge);
   aligner.setIdentityCutoff(param.ic);
   aligner.run();
   time(&timer2);
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << " for writting\n";
      return 1;
   }
   ouf << seq1.first << " x " << seq2.first << endl;
   aligner.displayAlignment(ouf);

   double seconds = difftime(timer2, timer1);
   cout << "Alignment took " << seconds << " seconds "
      << " and writte to file: " << outfile << "\n";

   return 0;
}

int selfAlign(const string &infile, const Progparam &param) {
   vector<DNA> allseq=slurpDNA(infile);
   cout << allseq.size() << " sequences for pairwise alignment\n";
   string outfile="test_pairwise.txt";
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << " for writting\n";
      return 1;
   }
   time_t timer1, timer2;
   time(&timer1);
   int backwardCnt=0;
   for (int i=0; i<allseq.size()-1; ++i) {
      //cerr << "working on i=" << i << endl;
      for (int j=i+1; j<allseq.size(); ++j) {
         ouf << "i,j=" << i << "," << j << ' ';
         Braboualn aligner(allseq[i].toString(), allseq[j].toString(), 
                  param.m, param.s, param.go, param.ge);
         aligner.setIdentityCutoff(param.ic);
         aligner.run();
         ouf << allseq[i].getName() << " x " << allseq[j].getName() << endl;
         aligner.displayAlignment(ouf);
         double forwardIdentity=aligner.getIdentity();
         if (forwardIdentity < 0.5) {
            ++backwardCnt;
            ouf << endl << string(20, '=') 
               << " Reverse Complement " << string(20, '=') << endl;
            DNA dnarc=allseq[j].revcompCopy();
            Braboualn rcaligner(allseq[i].toString(), dnarc.toString(), 
                        param.m, param.s, param.go, param.ge);
            rcaligner.setIdentityCutoff(param.ic);
            rcaligner.run();
            ouf << allseq[i].getName() << " x " << dnarc.getName() << endl;
            rcaligner.displayAlignment(ouf);
            if (rcaligner.getIdentity() > forwardIdentity) {
               allseq[j]=dnarc;
            }
         }
      }
      cout << endl;
   }
   time(&timer2);
   double seconds = difftime(timer2, timer1);
   cout << allseq.size()*(allseq.size()-1)/2 << " forward and "
      << backwardCnt << " alignments took " << seconds << " seconds "
      << " and writte to file: " << outfile << "\n";
   return 0;
}

string extractSeqid(const string &txt) {
   string id=txt.substr(1);
   size_t i = id.find("\t ");
   if (i != string::npos) {
      id=id.substr(0,i);
   }
   return id;
}

vector<DNA> slurpDNA(const string &fname) {
   ifstream inf(fname.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file: " << fname << endl;
      exit(1);
   }
   vector<DNA> result;
   DNA dna;
   string seq, line;
   while (dna.read(inf)) {
      result.push_back(dna);
   }
   inf.close();
   return result;
}

// more than one sequence in the file
// return a map of name->sequence
vector<pair<string,string> > readMultipleSequences(const string &fname) {
   ifstream inf(fname.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file: " << fname << endl;
      exit(1);
   }
   vector<pair<string,string> > result;
   string header, seq, line;
   getline(inf, line);
   while (!inf.eof() && line[0] == '>') {
      header=line;
      getline(inf, seq);
      getline(inf, line);
      while (!inf.eof() && line[0] != '>') {
         seq += line;
         getline(inf, line);
      }
      result.push_back(make_pair(extractSeqid(header),seq));
      //cout << "#" << result.size() << endl;
      //cout << header << " " << seq << endl;
   }
   inf.close();
   return result;
}


pair<string,string> readSequence(const string &fname) {
   ifstream inf(fname.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file: " << fname << endl;
      exit(1);
   }
   string header, seq, line;
   getline(inf, header);
   getline(inf, line);
   while (!inf.eof()) {
      seq += line;
      getline(inf, line);
   }
   inf.close();
   return make_pair(header.substr(1), std::move(seq));
}

string makeOutFileName(const string &file1, const string &file2) {
   string::size_type i = file1.rfind('.');
   string tmp;
   if (i == string::npos) tmp=file1;
   else tmp = file1.substr(0,i);
   tmp += "_";
   i = file2.rfind('.');
   if (i == string::npos) tmp += file2;
   else tmp += file2.substr(0,i);
   tmp += ".aln";
   return tmp;
}
