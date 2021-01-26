#include <iostream>
#include <gtest/gtest.h>
//#include <fastq.h>
//#include "umicooker.h"
//#include "helper.h"
#include <dynalnt.h>
//#include <api/BamReader.h>
//#include <api/BamAlignment.h>
#include "../minhash.h"
#include <iterator>
#include <string>
#include <kmerbase.h>

using namespace std;
using namespace orpara;
//using namespace BamTools;

#define KVALUE 19
#define SVALUE 25

// to be finished
class MinhashTest : public testing::Test {
   protected:
      MinhashTest() : bamf("/prednet/data03/OutputByRun03/pipeline_test/180816_NB501494_0092_AH3KHKBGX7/lbwfresult/cfDNA_PreValidation_inter_20180816/test/X15ng_60per_Rep1_merged_svmaterial.bam")
      { }
      virtual ~MinhashTest() { }
      virtual void SetUp() { }
      virtual void TearDown() { }
      const string bamf;
};

TEST_F(MinhashTest, repeatseq) {
   //vector<uint64_t> sketch1{0,2,8,32,128,515,2060,8240,32960,131840,527360,2109440,8437760,33751040,135004160,540016640,2160066560,3221225472,8640266240};
   //vector<uint64_t> sketch2{0,2,8,32,128,515,2060,8240,32960,131840,527360,2109440,8437760,33751040,135004160,540016640,2160066560,8640266240,34561064960};
   string seq1("GGGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGGGGAAAAAAAGGAAAAAAAAAAATAAAAATATATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAATAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA");
   string seq2("ACGGAGGTGACCTGAGTCCTGAAGGCGGAGGTTGCAGTGAGCCAAGATGGCACCACTGCACTCCAGACTGGGAGACAGAGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAATAAAAAAAAAAAAAATTTAA");
   Sketch<KVALUE, 60> sk1(seq1);
   cerr << "sketch for " << seq1 << endl;
   sk1.show(cerr) << endl;
   Sketch<KVALUE, 60> sk2(seq2);
   cerr << "sketch for " << seq2 << endl;
   sk2.show(cerr) << endl;
   float jsimilarity = sk1.jsim(sk2);
   cerr << "jsim=" << jsimilarity << endl;
   ASSERT_TRUE(jsimilarity < 0.25);
}

TEST_F(MinhashTest, falsepositive) {
   /*
sketch by C function Jcard similarity: 0.66667
similarity=0.6129
TTTTAACACCAACACCCCCTACCACACCCCTTTTACTGTCTGGGTTTTATTCCTGTGGTGTTCGGATATAAAAATTTTAAAAAAAAAAAAAAACAAAAACAAACAAACACACACACACACACACACAACCTTCTTTTTTTTTTATTTTGAG MBC=TCTTACGG+CCGTAAGA
AGCTGGGACTACAGGTGTGAACCATTATGCCTGGCTGGAAAGCTACTTTAGAAATGTTTTTTTTTAAAAAAAAAAAACAAAAACAAACAAACACACACACACACACACACACACACACACACACACATACCAAAATATTTATTTATGCATT MBC=TCGTCTGA+CGGTCATA
seq1 x seq2  46-144/151 | 24-145/151
Score=1054 gap length: 25 2 num gaps: 4 2 idencnt=80 simcnt=0 alnlen=124 identity=0.64516 similarity=0.64516
47        11        65        72        79        89        99        109
+         +         +         +         +         +         +         +
TTATTCCTG--TGGTGTTCGGA---TATAAAAAT---TTTAAAAAAAAAAAAAAACAAAAACAAACAAACACACACACAC
|||| ||||  |||    |  |   || |||  |   |||  |||||||||||| |||||||||||||||||||||||||
TTATGCCTGGCTGGAAAGCT-ACTTTAGAAATGTTTTTTTTTAAAAAAAAAAAA-CAAAAACAAACAAACACACACACAC
+         +         +         +         +         +         +         +
25        35        21        54        64        74        83        93

119       91        101       132       142       13
+         +         +         +         +         +         +         +
ACACACACA-----------------ACCTTCTTTTTTTTTTAT
|||||||||                 |||    | ||| |||||
ACACACACACACACACACACACACATACCAAAATATTTATTTAT
+         +         +         +         +         +         +         +
103       113       123       133       143       13
*/
   string seq1("TTTTAACACCAACACCCCCTACCACACCCCTTTTACTGTCTGGGTTTTATTCCTGTGGTGTTCGGATATAAAAATTTTAAAAAAAAAAAAAAACAAAAACAAACAAACACACACACACACACA    CACAACCTTCTTTTTTTTTTATTTTGAG");
   string seq2("AGCTGGGACTACAGGTGTGAACCATTATGCCTGGCTGGAAAGCTACTTTAGAAATGTTTTTTTTTAAAAAAAAAAAACAAAAACAAACAAACACACACACACACACACACACACACACACACA    CACATACCAAAATATTTATTTATGCATT");
   Sketch<KVALUE, 60> sk1(seq1);
   cerr << "sketch for " << seq1 << endl;
   sk1.show(cerr) << endl;
   Sketch<KVALUE, 60> sk2(seq2);
   cerr << "sketch for " << seq2 << endl;
   sk2.show(cerr) << endl;
   float jsimilarity = sk1.jsim(sk2);
   cerr << "jsim=" << jsimilarity << endl;

   //ASSERT_TRUE(jsimilarity < 0.25);
   Sketch<KVALUE, 60> sk3(seq1, KmerBase<KVALUE>::base2intC);
   Sketch<KVALUE, 60> sk4(seq2, KmerBase<KVALUE>::base2intC);
   jsimilarity=sk3.jsim(sk4);
   sk3.show(cerr) << endl;
   sk4.show(cerr) << endl << "jsimC=" << jsimilarity << endl;

   Sketch<KVALUE, 60> sk5(seq1, KmerBase<KVALUE>::base2intT);
   Sketch<KVALUE, 60> sk6(seq2, KmerBase<KVALUE>::base2intT);
   jsimilarity=sk5.jsim(sk6);
   sk5.show(cerr) << endl;
   sk6.show(cerr) << endl << "jsimT=" << jsimilarity << endl;
}

/*
TEST_F(MinhashTest, baminput) {
   //return;
   BamReader reader;
   reader.Open(bamf);
   BamAlignment* ba = reader.next();
   map<Sketch<KVALUE, SVALUE>, BamAlignment*> store;
   map<Sketch<KVALUE, SVALUE>,BamAlignment*>::iterator it, itt;
   string mbc1, mbc2;
   int cnt=0;
   int every=5000;
   while (ba != nullptr) {
      if (cnt % every == 0) {
         cerr << " reading #" << cnt << endl;
         every += 5000;
      }
      Sketch<KVALUE, SVALUE> mh(ba->getQuerySequence());
      if ((it=store.find(mh)) != store.end()) {
         it->second->GetTag("BC", mbc1);
         ba->GetTag("BC", mbc2);
         //it->first.show(cout) << " ";
         //mh.show(cout) << endl;
         //cout << "Almost identical sequence found before:\n"
         //   << mbc1 << "  " << it->second->getQuerySequence() << endl
         //   << mbc2 << "  " << ba->getQuerySequence() << endl;
         delete ba;
      }
      else {
         store.insert(make_pair(std::move(mh), ba));
      }
      ++cnt;
      //if (cnt > 2000000) {
      if (cnt > 20000) {
         cerr << " stop testing\n";
         break;
      }
      ba = reader.next();
   }
   reader.Close();
   cerr << store.size() << " unique sequences trying to find similar sequences\n";
   SimpleScoreMethod sm(20, -20, -50, -2);
   Dynaln<SimpleScoreMethod> dyaln(sm);
   for (it=store.begin(); it != std::prev(store.end()); ++it) {
      for (itt = std::next(it); itt != store.end(); ++itt) {
         float js = it->first.jsim(itt->first);
         //cout << "js=" << js << endl;
         if (js > 0.61) { // very similar sequences
            if (!it->second->GetTag("BC", mbc1)) {
               cerr << "Failed to get BC tag from Bam\n";
            }
            if (!itt->second->GetTag("BC", mbc2)) {
               cerr << "Failed to get BC tag from Bam\n";
            }
            if (mbc1 == mbc2) {
               cout << "Same MBC\n";
            }
            cout << "similarity=" << js << endl;
            //it->first.show(cout) << endl;
            //itt->first.show(cout) << endl;
            cout << it->second->getQuerySequence() << " MBC=" << mbc1 << endl
               << itt->second->getQuerySequence() << " MBC=" << mbc2 << endl;
            DNA dna1("seq1", it->second->getQuerySequence());
            DNA dna2("seq2", itt->second->getQuerySequence());
            dyaln.setseq(dna1, dna2);
            dyaln.runlocal();
            dyaln.printAlign(cout, 80);
            cout << string(60, '-') << endl;
            Sketch<KVALUE,SVALUE> sk1(it->second->getQuerySequence(), KmerBase<KVALUE>::base2intC);
            Sketch<KVALUE,SVALUE> sk2(itt->second->getQuerySequence(), KmerBase<KVALUE>::base2intC);
            float newjs = sk1.jsim(sk2);
            cout << "sketch by C function Jcard similarity: " << newjs << endl;
         }
      }
      delete it->second;
   }
   delete it->second;
}
*/


