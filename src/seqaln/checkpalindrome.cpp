#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <bitset>
#include <cstring>

#include "bioseq.h"
#include "kmert.h"
#include "fastq.h"
#include "kmer.h"

using namespace std;
using namespace orpara;

int checkOne(string name,string title,int len,const string &seq, vector<double> &freq,vector<double> &rcfreq) {
   Kmert<7> kmer6(seq); // actually 6
   pair<int,int> loopReg;
   int cutPoint = -1;
  // cout<<seq<<endl;
  
   if (kmer6.isPalindrome(loopReg, freq, rcfreq)) {

     cout << name << " " << title  
            << " length: " << len<< endl;



     

      /*vector<int> dist2end;
      //cout << "loop region: " << loopReg.first << "-" << loopReg.second << endl;
      dist2end.push_back(loopReg.first);
      dist2end.push_back(loopReg.second);
      dist2end.push_back(seq.length() - loopReg.first);
      dist2end.push_back(seq.length() - loopReg.second);
      int mind=99999999;
      int mini=0;
      for (unsigned int i=0; i<4; ++i) {
         if (abs(dist2end[i]-1500) < mind) {
            mind=abs(dist2end[i]-1500);
            mini = i;
         }
      }
      //cout << "mind: " << mind << " mini: " << mini << endl;
      if (mini == 1 || mini == 3) {
         cutPoint=loopReg.first;
      }
      else {
         cutPoint=loopReg.second;
      }*/
   }
   //else {
   //   cout << "is not palindrome\n";
   //}
   return cutPoint;
}

string nameOutfile(const string &infname) {
   string::size_type i = infname.rfind('.');
   string outfname = infname;
   string fext = "fas";
   if (i != string::npos) {
      outfname = infname.substr(0, i);
      fext = infname.substr(i+1);
   }
   outfname += "_nopalindrome." + fext;
   return outfname;
}

void usage() {
   cerr << "checkpalindrome -o outfile_name infput_name\n"
      << "input_file must be fastq or fasta file\n";
   exit(1);
}

void checkFasta(ifstream &inf, ofstream &ouf, vector<double> &freq, vector<double> &rcfreq) {
   DNA dna;
   int cp;
   while (dna.read(inf)) {
      if ((cp=checkOne(dna.getName(),dna.getTitle(),dna.length(),dna.toString(),freq,rcfreq)) != -1) {
         cout << dna.getName() << " " << dna.getTitle()  
            << " length: " << dna.length() << endl;
         cout << "cutting sequence at " << cp << endl;
         DNA dna1=dna.subsequenceWithName(0, cp);
         DNA dna2=dna.subsequenceWithName(cp);
         ouf << dna1 << dna2;
      }
      else {

        // cout << dna.getName() << " " << dna.getTitle()  
         //   << " length: " << dna.length() << endl;

         ouf << dna;
      }
   }
}

void checkFastq(ifstream &inf, ofstream &ouf, vector<double> &freq, vector<double> &rcfreq) {
   Fastq dna;
   int cp;
   int numcut=0;
   int total=0;
   while (dna.read(inf)) {
      ++total;
      if ((cp=checkOne(dna.getName(),dna.getTitle(),dna.length(),dna.getSequence(),freq,rcfreq)) != -1) {
         cout << dna.getName() << " " << dna.getDescription()  
            << " length: " << dna.length() << endl;
         cout << "cutting sequence at " << cp << endl;
         Fastq dna1=dna.sub(0, (unsigned)cp);
         Fastq dna2=dna.sub((unsigned)cp, dna.length()-1);
         ouf << dna1 << dna2;
         ++numcut;
      }
      else {

        // cout << dna.getName() << " " << dna.getDescription()  
          //  << " length: " << dna.length() << endl;
         ouf << dna;
      }
   }
   cout << numcut << " palindrome cut out of " << total << endl;
}

bool isFastq(const string &fn) {
   if (fn.rfind('.') != string::npos) {
      if (fn.substr(fn.rfind('.')+1) == "fastq")
         return true;
   }
   return false;
}

/**
 * need a lots more work to do to make it
 * a useful program. Right now is it still 
 * in development stage
 */
int main(int argc, char* argv[]) {
   //string infile="palindrome.fas"; // a few test sequences that are palindrome
   //string outfile="nopalindrom.fas";
   string infile, outfile;
   bool isFasta=true;
   int i=1;

 string libfile="silva123.fas";

      ifstream inf1(libfile);

      if (inf1.fail()) {

         throw runtime_error("failed to open reference library file: " + libfile);

      }

      DNA dna;
     
      int c=0;

      KmerCount counter(5);

   	while (dna.read(inf1)) {

      		counter(dna.getSequence(),c);

        }
    
   vector<double> ref_freq=counter.getFrequency();

   vector<double> ref_rcfreq=counter.getRCFrequency();



   for(int i=0;i<ref_freq.size();++i){

      ref_freq[i]=0;



   }

    while (i < argc) {
      if (!strcmp(argv[i], "-o")) {
         outfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-q")) {
         isFasta=false;
      }
      else {
         infile = argv[i];
      }
      ++i;
   }

   if (infile.empty()) usage();
   if (outfile.empty()) outfile = nameOutfile(infile);

   ifstream inf(infile);
   if (inf.fail()) {
      cerr << "failed to open " << infile << endl;
      return 1;
   }
   ofstream ouf(outfile);

   if (isFastq(infile)) {
      isFasta = false;
   }

   if (isFasta) {
      checkFasta(inf, ouf,ref_freq,ref_rcfreq);
   }
   else {
      checkFastq(inf, ouf, ref_freq,ref_rcfreq);
   }


   return 0;
}

