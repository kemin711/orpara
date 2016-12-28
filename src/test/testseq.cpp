#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "codon.h"
#include "bioseq.h"
#include "fastq.h"

using namespace std;

void testcodon() {
   cout << "Test codon table ...\n";
   codon::setCodonFile("/home/zhouke/etc/codontable.txt");
   DNA::setCodonTable(12);
}

void testFindORF(const string &file) {
   cout << "testing findORF from file: " << file << "\n";
   if (file.empty()) {
      cerr << "no input file\n";
      return;
   }
   bioseq seq;
   seq.read(file);
   vector<Range> orfs=findAllORFIndex(seq.toString(), 1, 70, 160, 120, 10);
   cout << orfs.size() << " ORFs found from RNA seq len=" << seq.length()
      << "\n";
   if (orfs.size() == 0) {
      cout << seq << endl;
   }
   for (unsigned int i=0; i<orfs.size(); ++i) {
      cout << orfs[i] << " " << orfs[i].length() << endl;
   }
   string pep; int bb, ee, ff;
   longestORFPlus(seq.toString(), pep, bb, ee, ff);
   cout << "RNA\n" << seq << endl
      <<"Pep:\n" << pep << endl
      << " " << bb << "-" << ee << " " << ff << endl;
   cout << "Tesing findORF() done\n";
}


void testLongestNMissingORF() {
   cout << "test LongestNMissingORF ...\n";
   string rna="GCAAACAAACAAAGCGGATCGAAAGACAAGGGACAGAACAAAAGACACATCCCCCTCCGGAGCCGTAAAAGGCAGGAAAATCTCCATTACCGCCGCTCCTTTCCCAGCGCGCTAGTTCAAGTGATAAACAGAGGCTAAGTTTCTATTAGCTTGCGACATTGCTGAACGTGAAGGGAACAAACGCACCTACAATGTCAAGGCTCAAAAATGCCCTTTTGAAACGCCAACTGGATGACATCGTCGAGCCAGTGGTTGCATCGATAAAGAGGCTTGGAAGGGTCCTTGCGTTTACCACCAAGAGTAACCTGAGCCCCTGGTGGGGCAACTGTATTACAGCACTCCGTTTCGAACCTTGTCAGGTTTTCAACCATGATTGTTCCGATCTTCTTTTTGGTGTGCGCCACGTAGGTGGCCGAGAGGAGAGGGTCCTTCGTGTCTCGCTTTTGGTAGGCCATGCCTTGGAGGATATTTCCCACGACGTGATACATTAGTCCGTTCCGAGTCTGTTCATCTGCGCTGGGCTGCGTTTCGATATAGATAGCAACATGGTTTGGCTGGCCAATGTATTCCGCAGTAAAGACTTCGTAAGACATAGTAGAAGGGTAATAGAAACTGGAAATCTATGAAGAGTATAGATGAAATTAAGTGAATAAGTTGGGAAGAACGAGAGTTTGAGTATATAGTTGAGTAGAATCAGGTCATTCTCGTGAGACTCGACATCCTTTAAAAGCTTGCCTTCCAGCAGATAGCAAGCTGTGGTAGTGAGGCGATGCATAGAAGATGGCGCCTACTACCTAGTCACCTACTTCCAACAGGATAGGGACCCATCGTGGCGGCTACCTTGTTGTGAAGCCCATGAATATCAAGCTGACCAATCTGGAGACTCGTCTAAGCACACACACTAATCCTGAAACACTTCTTGTCCCAGAGAATTAAATAGAAAAAATGAGCTATCATGCTAATTGGTTTAGTATATGAATGGTGTATGCTTGTAAGCGCCACCTCTAATAATAGCTGCAACAAGGAGCTCGTTTGCAGTAGGACTCCAAGTAGCCGCTGGAAATGGGCCCGACTAGCTCGATAAAACCTGTAACCATTGAAGCTTCCACTTGAACTTTGGGGAATTGATTGGCTCTCGTGGGCAATGCTTGATAGCTTGGAGTGTGTGGCCACCTGTTGATGTCCATGCCTGGGATCCCTTGGCATGCTTAGTAGGCCTATATGTCTCCAAATCGTGACAATCTGATGCATTGGGGGGTGGAATGAGTGGTTCTGCTGGTCGGTTGGAGCCAAGCGAGGCGCCAAAGGCCTCGAGGGTGCTGCGAGCAAGCTATGCCGCAAAAAATGACTGCGTAATCGATGGCCCAACTCTAAGCGCCGGCTGAATAGTTCCCGTGGCAGATCCGCCAATTGACTTTCCAATCAGGACCCTCTGGGTTTTTGAAATCGTTCCCATGGACGTAAATTACATCAGGATCTCGTTTGTTCTGTGCCTCGGCGCAGGCCGTTAATGACGGATATTGCACACGCGTACTCCGAAGCAATTTTTGAATAGCTGCCCGAGGCTCTAGGTCCACTACAGGCTCAAGTATGTCATCTAAGTATATAGAGGATCGGAGAGTTCGGGGAGGACGACCCGAGCGGTAGGACTGAGTGTGTAGGCCTTGCGGCCACACCTTGCGCTCCGGAGCCACATGGTCGAGAGATCGGGGCAGTACGCCAGCACGGGCTGGAAATCACCGACAGACATCCACGCCAGGAGCAAGACCCTTATGGTGGTGTATCTATGTTCTTGCGGCCCCATTGGTGTGCAAACCTTTTAACTCTTCTAGATCCAATTTTCGTGGTCCTCCCATAGAAATTTTTGGATACCAGAATTTTCTGGCCTGGTATTACTTACTTCTCGGTTATACTCTCATACTAGAAAGGGTTGTGCTGGGAAGAGCTTGATCCT";
   pair<int,int> Nterm1, full1, Nterm2, full2;
   string nspep1, pep1, nspep2, pep2;
   longestNoStartORFPlus(rna, Nterm1, nspep1, full1, pep1);
   cout << "no start pep: " << Nterm1.first << "->" << Nterm1.second 
      << endl << nspep1 << endl << "full pep: " << full1.first 
      << "->" << full1.second << endl << pep1 << endl;
   pair<int,int> best;
   string pep;
   if (full1.second-full1.first > Nterm1.second-Nterm1.first) {
      cout << "use full length ORF\n";
      best=full1;
      pep=pep1;
   }
   else {
      pep=nspep1;
      int frame=Nterm1.first;
      cout << "Frame: " << frame << endl;
      best=Nterm1;
      best.first=0;
   }
   string rcrna=rna;

   reverseComplementInPlace(rcrna);
   longestNoStartORFPlus(rcrna, Nterm2, nspep2, full2, pep2);
   cout << "On reverse strand\nno start pep: " << Nterm2.first 
      << "->" << Nterm2.second << endl << nspep2 << endl
      << "full pep: " << full2.first << "->" << full2.second << endl
      << pep2 << endl;
   
   longestNoStopORFPlus(rna, Nterm1, nspep1, full1, pep1);
   cout << "no stop ORF:\n" << Nterm1.first << "->" << Nterm1.second
      << endl << full1.first << "->" << full1.second << endl
      << pep1 << endl;
   string shortrna="CAAACCACCCTCTTCATCGCGTTGTCAGAGCGATTCCCCGAAGAGCGCTAGGAACTCGCGACAATGCTTATCCCCAAGGAAGACCGCAAGAAGATCCACGAGTACCTCTTCCGCGGTATGTGAACACCCATATATCCACCGCCCCGTCCCCGGGTCGCGCAACCTTCGATCGATCCATCCGGTATCACTATGGAGGAGGAACAGTCTCGCTGACGGTTAGAACAGAGGGTGTGCTCGTGGCCAAGAAGGACTTCAACCTTCCCAAGCATGGCGACATTGACACCAAGAACCTCTACGTGATCAAGGCCCTGCAGTCCCTGACCTCCCGCGGTTATGTCAAGACCCAGTTCTCGTGGCAGTACTACTACTACACCCTCACCCCCGAGGGTCTTGACTACCTCCGTGAGTGGCTCCACCTCCCCGCTGAGGTTGTCCCCGCCACCCACATCAAGCAGCAGCGTTCCCACGCTCCCCCCCGTGGCATGATGGGTGGTGAGGACCGTGAGCGTCGTCCCCGTGCTCCTCGTGAGGGTGGCTACCGCCGCCGCGAGCAGGAGAACAAGGAGGGCGGTGCCCCCGGCGAGTTCGCTCCCAGCTTCCGTGGTGGATTCGGCCGTGGCCGTGGTGCTCCCTCCTCCTAAGGGGTGTCGCCTGCTGGTCTTCCAAG";
   longestNoStartORFPlus(shortrna, Nterm1, nspep1, full1, pep1);
   cerr << "Test longestNoStartORFPlus\n";
   cerr << Nterm1.first << " " << Nterm1.second << endl << nspep1 << endl;
   cout << "test LongestNMissingORF done!n";
}

void testLongestPlus() {
   cout << "testing LongestPlus() ...\n";
   string rrr="CTTGGTCGGTCTTGCTGCTCTCAAGGACTTCTTCCAACCGTGCTGACATAAGGGCCTCTTGCTCTACCGAGGCCGAAACGGCTGTGCACCCCGCCCTGGAACGTGCCCCTTCTTATATAACCGATCCATAGTAGTCGCATGTCTGCGTCTTCCTCCAGTATCTCCTCTCGCATTGTAAGAGCATCAATCGCCCGATCCTCTTTCCTCACGTGACCAACTGCTGCAGCTTTGATCGCCGTCCTAAAGCATTCAGTGAGTCACCGATCGTTGTGAGATACACCGGCCCATAGTATATTGCATTCAGACCCTACCTGTCTTCCCCCACAGCGTAAATAAGATTCCTGTGAGTTGATGTGTGTGCTGTGCTCTGAATATTCCTGTATCCCTTGCTTGCTCTTGACCACGCAGCAAAGTGCGCTTGTCCAAGGGACCCGGCGCAGTGTAAATAGCTGACACTATATGAAGCGTGCCAGCTGCAGATATATTGCATTGATACGCTGTTTATACCCCGCAAGGCATGGCCACGCAGCTCGTTTCCCTGGCGGAGGTGGAGAGGCTAAGTGCTTCGGTCGTGCGGATTTTGGGAGGAAACCCCGGCAAGGCTTCTTTCTTTTAGTTCACGCTACAAGG";
   string ppp;
   int bb, ee;
   longestORFPlus(rrr, ppp, bb,ee);
   cerr << "Result of longestORFPlus:\n" << ppp << endl << "orf index " << bb << " " << ee << endl;

   pair<int,int> Rfull, Rnoend;
   string Pfull, Pnoend;
   longestNoStopORFPlus(rrr, Rnoend, Pnoend, Rfull, Pfull);
   cerr << "result of longestNoStopORFPlus:\n"
      << "partial: " << Rnoend.first << " " << Rnoend.second << endl
      << Pnoend << endl
      << "full: " << Rfull.first << " " << Rfull.second << endl
      << Pfull
      << endl << "RNA length: " << rrr.length() << endl;
   cout << "testing LongestPlus() done!\n";
}

void testLongestORF() {
   cout << "testing longestORF() ...\n";
   string rnaseq="GCACAGCGGCGGGCGGTGCACGTGTGCTGCACTTCCCGCAGCCCGTCCGCATCTTCAGCGGCTTCAACAACCACGCCACCTGGGACCGATTTGACGAGCTCATGCAGCGCCACACCACACACTGGTGCTGCCGCAGCCCGCCAGACATGAAGGCCTACAACCTGACGGAGCGGCTGCAGCTGGTGGCGCTGCCGCCCGAGCGCTACAAGAGCCTGCCGCCGCTAGAGGCGCGAACCAGCTACCTGCACACACTGGGGCCATGGCAAGGGCTCAGCAAAGGGACGTAGCGCTGCCCACTGACAAGTGTGCGGTGGTTTCAAGAATACACTAGCGTATGGACAGCGTGGAAGATGCGATGATGGGAAGCAAAGGGCAAGTCCGAGAGTCGCGCAATGTTGCAGGAGTATGTGTGTAGATGGGTAGGTGCGTGTGTCTGGCGGTGCGGGGCTAGAGAAGCGAGAAGACAGTGTAGGACTTCATGGTGTGGGACTTGAACTAGGAAGGAATCTGACGGTGCTTAATGACTGGTGGATACCGGACGTCTTTTGGTGGCTCGCGGGACCAACCCGGTGGGCCGCCCAGACACGGAGGGGGAAAGGACACGCGAGCGCGGATGATGGGCCCCACACCCACTTCATCACCCCCACCCGCTGCTGAGCTGGCTACATACACTGCTTTTAAGTGTCATACTGGCCCTTGGCTAGTAAGCGACGGAAACCACGGGCTTCTTTGCTCTGTTATTGCCACACATCCGCGACGCTCGGCGTGGAGCTGCACACAGCGACCAGCCATTGCAACGATTTGGGCAGCTTCGGTCGCCAGGAGCTTGAGCGCGTGGAGACTTCTCTGGACGCGTGAGCAGCCTGCTTTCCACTTGATACTCCTTGCACCACCGAGCACCGACACCATGCAGCACACGCTCATGCGCGCTAATTCCTACCCAGCGGGCGCTGCCTGCAGCGGTCGCGCCTCGCAGCGCCGCGCGCAGCCCGTCACAGTGTGCAGCGGCGCCAGCCGCACGCCGCAGCAGGGCATCGCCCCGATGCGCGGCAGCTCGGCCTGCGGCGCCTCCGGCCGCCTGCCCGCGCCCTTCTCGTCCGCCGCCGCCTCGCGCCCCACCGCCGCTGGCCGCCGCGGCGCTGTGCGCGTGCAGGCCAACTGGGGCGCCCCGGTGGAGTTCCAGCCTGCAAAGGTTGTGTCCAACTCGCCCGCCGCCGCCGGCCCGCTGCACAAGGTGGTGATCGACGTGGGTGCGCCGCTGGCCGCCGGCTACACGGTGCCGGGCCAGTTCGTGCAGGTCAAGGTGGGGGACAGCAAGCCCGGCTT";
   Protein pep;
   int b,e;
   DNA rna(rnaseq);
   rna.longestORFForward(pep,b,e);
   cerr << b << " " << e << endl << pep << endl;
   cout << "testing longestORF() done!\n";
}

void testTranslate() {
   cout << "testing translate() ...\n";
   DNA dna("ACGTAAAATTTAAACCCTTTTAAAACACACGGAATCGCACGCTCACTAACATCGCGCACATTACGCATTTGCTCGCGCAATGTCATGCAGGTGCACATCATCAAGGCGGTGTACGAGTTCGGCGTGAACGTGATCCACTCTGACACGGACGTGGTCTGGTTCGGCGACCCGCTGCCATTCTTCCAACAGCAACTGCAGGAGAGACACGGAGGCGACGGCGGCGATGACGGCATCGGCCACAAAGTTCGGGCGAAGTGGGACGCAGCAGCGTCAGAGGTTGGGGCTGCTGGGGCCACCCCCGCGGCTGCGGACGCGGCTGCGGCAGGGTCGTCGGGCTGGTGGGCGCCGCCGCATGTGGCGGTGGCCACCGACTCGGTGTCCACTATGAACCGCCGCGGCGACCGCGGCCTGGAGGACAGCCCGCAGCCGTACGCGCCCATCAACACAGGTAGGTAGGCGCATCATGGCCAGGGCGACAACCAACAACGCGAGGGTGCATGTTTGGACTGCTAGCTACCGTAGTTCAAGGGTGCGGGCAGGTGGGCGGCCCCGCCCACAGCGCGGCTCGGGTGCGGACTGACCCGACTGGCCCCAATGCATGCATGGCAAGATCCATGCATAGCAACCAGCCACCACCATGTCAGCACCCAACCGCCAGGCCGCCACACCTGCACATGCCAACACGCGTGACACGACCCATGCATTGGCCCGCCCGGTTAGCTCTCCCCCTTCCACAGCTCCATCCACGCATCCTCAACCCTATGCACTACACTACACTGCCTGCAGGCATCTATTTCCTGCGGCAGTGGTCGGGCGGGGCGGCATTCCTGGATGAGTGGCTGTCGTGGCAGGAGCGGGACGTGGGTCACGACCAGGTAGGTGCTGATGCGATTAATGAGGCCGCGGGCGGGAGGGGGGCCGTTCGGAGGACCGAATGATGCCTTCTTTGGGAGGGGCCGGACGCGAGTGTCGTGCATGGGGGCGGGTGAGGTCAACTGCTGCTGGGCCTTAGGACTACCGCACCCACACACATTTATGCGCCCTCACCTGTATGTGGGCCTGCCGACACGTGCGTGCACCCACGTACATGCACACACGTGCCCACATATGCACAGGACGGTTTCAACACACTGGCCCGGGGCTTCTTCTTCCACCGCGACCCCGACCTGCGCCTGCCGGTGTTCCCGGATGCACACACTCGCATAACACTGCCCCGGGGCAGTACCAGCGGCAGCGGCACTAATGATGGCAGTGAACGGGACGTGGGCGGCGTTGGTGACGGCAACACGACGGGCGAGTGGGCGGCACCGCGCACCTTCCGTGCCGCCTACAGCAACACCACGGGCGTGTCCTTCCTGCCCGCATCCATGTTCGGCAACACCTACACCTACGTCAATGCGCGGCTGTGGGAGGTGGGAACAGGGACAGGGGGGGGTGCGGGGCTGGGAGGGAGACCGCAGCGGGAGAGGCTGCGGACTTGCAGGGCGGGTGGTCTGCGGCGGTCTTCACATGCATACGGTACACAGCCGACAGGTGTACCATCCTGACATGTATCATGCCCCATGCCTCATGACCCCACCTCTTCGCCACCTATCATGTGTTTCGCAGAAGCTAAAACACCCGCTGTATGCCATCCACTGGGTGTGGGGTGGGAGCACACTTGAGAGCAAGCGCCAAAATATGCGCGATGCCATGAAGTTCCACGATGAGCCCGACTACTACACCTCGCCAAACCTGGTCACCTTTGACCTGGACCTGCTGCCGGTGCGTGCATGCCGTACGGTTTGTTGTGCCATGCGAATGATGCGACCAATTCACCCCCATTCACATAATTCACGCCGCCCCTACGTACCGCCGCGTGCGCCCCTACTGCTCTGGCGCCCATGGCTCGATGCAGGTGCCGCCCACCTTCAACTCCTGGTTCAGTACCGAGCACATGATCCGCTTCCACGTCCAAGCCGCCAATTATCAGCTGCAGCAGGCTTACTACGCCTTTGCCATCGCGCTCATTGCCAACCGCACGCTGGTCATGCCGCGGGTACGTACGCGTGCGGT");
   int begin, end;
   Protein pp = dna.longestORFForward(begin, end);
   //cout << begin << "-" << end << endl;
   //cout << "From reverse strand\n";
   Protein pp2 = dna.longestORFReverse(begin,end);
   cout << begin << "-" << end << endl;
   Protein ll=dna.longestORF(begin, end);
   cout << "longest EST from both directions\n";
   cout << ll << endl
      << begin << " " << end << endl;
   cout << "testing translate() done!\n";
}
// 1  2  3  4  5  6  7  8  9  10 11 12 13 14                               
// TGCATTTTAAATTTGGGAAAATTTTGTGTGCCTTAGCGTGC
// ACGTAAAATTTAAACCCTTTTAAAACACACGGAATCGCACG
// 0    5    10        20        30        40
// L=41

void testLoad(const string &file) {
   cout << "test load() ... \n";
   map<string,string> store;
   loadFastaIntoMap(file,store);
   cout << "first sequence\n";
   cout << store.begin()->first << endl;
   printFasta(cout, store.begin()->second.substr(0, 2000));
   cout << endl;
   cout << "test load() done! \n";
}

void testDNAQual() {
   string qfile="test.fastq";
   ifstream inf(qfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << qfile << endl;
      exit(1);
   }
   vector<DNAQual> reads;
   Fastq fastq;
   // read the objects from file
   while (fastq.read(inf)) {
      reads.push_back(DNAQual(fastq.getName().substr(1), fastq.getSequence(), fastq.getQuality()));
      if (fastq.hasDescription()) {
         reads.back().setTitle(fastq.getDescription());
      }
   }
   // save for later testing.
   DNAQualCount dqc(reads.back(), 99);
   DNAQualCount dqc2(reads.front(), 88);

   cout << "\ntesting getRevcomp() ...\n";
   const DNAQual *revcpl;
   for (size_t i=1; i<reads.size()-1; ++i) {
      cout << reads[i] << endl;
      cout << string(70, '-') << endl
         << " after calling getRevcom()\n";
      revcpl = reads[i].getRevcomp();
      cout << *revcpl << endl;
      cout << string(70, '*') << endl;
   }
   cout << string(80, '=') << endl;
   cout << "\nRound 2 testing getRevcomp()\n";
   for (size_t i=0; i<reads.size(); ++i) {
      cout << *(reads[i].getRevcomp()) << endl;
   }
   cout << "testing revcomp() after calling getRevcomp()\n";
   // reverse complement inplace
   for (unsigned int i=0; i<reads.size(); ++i) {
      reads[i].revcomp();
      cout << reads[i] << endl;
      cout << string(80, '!') << endl;
   }
   // now play with the last object.
   cout << "DNAQualCount object\n" << dqc << endl;
   cout << "Make a copy through revcompCopy()\n";
   DNAQualCount Drc = dqc.revcompCopy();
   cout << Drc << endl;
   dqc.revcomp();
   cout << "after revcomp() in place operation\n";
   cout << dqc << endl;

   cout << "before revcom() of a fresh object\n";
   cout << dqc2 << endl;
   dqc2.revcomp();
   cout << "after revcom() of a fresh object\n";
   cout << dqc2 << endl;
   cout << "doing recvomp() second time\n";
   dqc2.revcomp();
   cout << dqc2 << endl;
   cout << "Number of identical reads " << dqc2.getCount() << endl;
   cout << "end of DNAQual and DNAQualCoutn testing\n";
}

void testSubseq() {
   string qfile="test.fastq";
   ifstream inf(qfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << qfile << endl;
      exit(1);
   }
   vector<DNAQual> reads;
   Fastq fastq;
   while (fastq.read(inf)) {
      reads.push_back(DNAQual(fastq.getName().substr(1), fastq.getSequence(), fastq.getQuality()));
      if (fastq.hasDescription()) {
         reads.back().setTitle(fastq.getDescription());
      }
   }
   // save for later testing.
   vector<DNAQualCount> xxy;
   DNAQualCount dqc(reads.back(), 99);
   DNAQualCount dqc2(reads.front(), 88);
   DNAQualCount sb=dqc.subseq(1, 100);
   cout << "original:\n" << dqc << endl;
   cout << sb << endl;
   xxy.push_back(sb);
   DNAQualCount sb2=dqc2.subseq(2,150);
   xxy.push_back(sb2);
   cout << sb2 << endl;
}

void testFastq() {
   string fastqFile="test.fastq";
   ifstream inf(fastqFile.c_str());
   if (!inf) {
      cerr << "Failed to open " << fastqFile << endl;
      exit(1);
   }
   Fastq fastq;
   while (fastq.read(inf)) {
      cout << fastq.getAverageQuality() <<  " "
         << 1-Fastq::q2p(fastq.getAverageQuality())<< endl;
   }
}


int main(int argc, char* argv[]) {
   cerr << "tesing this library ...\n";
   /*
   string str("ksksksksksdd very long line or short line");
   str = "assign a new value";
   vector<string> vstr;
   vstr.push_back("ssslsl");
   vstr.push_back("second");
   vstr.push_back("third");
   unsigned i;
   for (i=0; i<vstr.size(); ++i) {
      cout << vstr[i] << " ";
   }
   cout << endl;
   cout << "Done with simple test\n";

   codon ct;
   char cdo[3] = {'A', 'A', 'A'};
   char AA = ct[cdo];
   cout << "Amino acid " << AA << " has codon: AAA\n";

   testcodon();
   testFindORF("testseq.fas");
   testLongestNMissingORF();
   testLongestPlus();
   testLongestORF();
   testTranslate();
   testLoad("testloadmap.fas");
   cout << "reverse complement of A is " 
      << DNA::toRevcompBase('A') << endl;
   */
   testFastq();
   testDNAQual();
   testSubseq();
   return 0;
}

