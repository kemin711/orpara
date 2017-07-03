#include "bioseq.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "scorematrix.h"
#include <gtest/gtest.h>

using namespace std;
using namespace orpara;


class ScoreMethodTest : public testing::Test {
   protected:
      ScoreMethodTest() 
         : smbase(), ssm(20, -19), psm(), nsm()
      { }
      virtual void SetUp() { }
      virtual void TearDown() { }

      // members, 4 different matrix
      ScoreMethod smbase;
      SimpleScoreMethod ssm;
      ProteinScoreMethod psm;
      NucleicScoreMethod nsm;

};

TEST_F(ScoreMethodTest, defaultMatrix) {
   cout << "base scoremethod\n";
   int simval = smbase.lookup('A', 'A');
   ASSERT_EQ(simval, 10);
   cout << "A|A=" << simval << endl;
   simval = smbase.lookup('A', 'G');
   ASSERT_EQ(simval, -10);
   cout << "A|G=" << simval << endl;
   cout << "match: " << smbase.getMatch() 
      << " misMatch: " << smbase.getMismatch()
      << " gap open: " << smbase.getGapOpen()
      << " gap extend: " << smbase.getGapExtend() 
      << endl << endl;

   cout << "SimpleScoreMethod\n";
   simval = ssm.lookup('A', 'A');
   ASSERT_EQ(simval, 20);
   cout << "A|A=" << simval << endl;
   simval = ssm.lookup('A', 'G');
   ASSERT_EQ(simval, -19);
   cout << "A|G=" << simval << endl;
   cout << "match: " << ssm.getMatch() 
      << " misMatch: " << ssm.getMismatch()
      << " gap open: " << ssm.getGapOpen()
      << " gap extend: " << ssm.getGapExtend() 
      << endl << endl;

   cout << "ProteinScoreMethod\n";
   simval = psm.lookup('A', 'A');
   cout << "A|A=" << simval << endl;
   simval = psm.lookup('A', 'G');
   cout << "A|G=" << simval << endl;
   cout << "match: " << psm.getMatch() 
      << " misMatch: " << psm.getMismatch()
      << " gap open: " << psm.getGapOpen()
      << " gap extend: " << psm.getGapExtend() 
      << endl << endl;

   cout << "NucleicScoreMethod\n";
   simval = nsm.lookup('A', 'A');
   ASSERT_EQ(simval, 5);
   cout << "A|A=" << simval << endl;
   simval = nsm.lookup('A', 'G');
   ASSERT_EQ(simval, -4);
   cout << "A|G=" << simval << endl;
   cout << "match: " << nsm.getMatch() 
      << " misMatch: " << nsm.getMismatch()
      << " gap open: " << nsm.getGapOpen()
      << " gap extend: " << nsm.getGapExtend() 
      << endl << endl;
}

TEST_F(ScoreMethodTest, blosum) {
   ProteinScoreMethod pm("blosum62.50");
   cout << "testing protein matrix blosum62.50\n";
   ASSERT_EQ(pm.lookup('G', 'G'), 6);
   ASSERT_EQ(pm.lookup('W', 'W'), 11);
   cout << "G|G=" << pm.lookup('G', 'G')
      << " W|Y=" << pm.lookup('W', 'Y')
      << " W|W=" << pm.lookup('W', 'W')
      << " match: " << pm.getMatch() 
      << " misMatch: " << pm.getMismatch()
      << " gap open: " << pm.getGapOpen()
      << " gap extend: " << pm.getGapExtend() 
      << endl << endl;
}

// needs to remove text
// TODO: test failed, need to work this out!
/*
TEST_F(ScoreMethodTest, words) {
   ProteinScoreMethod pm;
   int neighborThreshold=10;
   float neighborFraction=0.75;
   vector<string> wds;
   pm.getWords(wds, 3);
   // too many to show just count
   //for (int i=0; i<wds.size(); i++) {
   //   cout << wds[i] << " ";
   //   if ((i+1)%10 == 0) cout << endl;
   //}
   cout << endl;
   cout << wds.size() << " words. vector<string> version tested fine\n";
   ASSERT_EQ(wds.size(), 32768);

   cout << "testing allwords() ...\n";
   char **aw=pm.allwords(2);
   ofstream outf("neighbors.txt");
   ofstream BAD("noneighbor.txt");
   //pm.showWords();
   map<int, int> stat;
   int numn;
   for (size_t i=0; i<pm.getNumberOfWords(); i++) {
      numn = pm.similarWord_debug(outf, aw[i], 2, neighborThreshold, neighborFraction);
      ++(stat[numn]);
      if (numn == 0) {
         BAD << aw[i] << endl;
      }
   }
   map<int,int>::const_iterator it;
   cout << "#number of neighbor | count at threshold=" << neighborThreshold
      << " fraction=" << neighborFraction << "\n";
   for (it=stat.begin(); it != stat.end(); it++) {
      cout << it->first << "\t" << it->second << endl;
   }
}
*/

TEST_F(ScoreMethodTest, nucleic) {
   NucleicScoreMethod nuc("NUC.4.4");
   cout << "min score: " << nuc.getMinScore() << endl;
   cout << "max score: " << nuc.getMaxScore() << endl
      << " gap open: " << nuc.getGapOpen() << " gap extend: " << nuc.getGapExtend()
      << endl;
   nuc.show();
   string base("ACGTN");
   for (size_t i=0; i<base.size(); ++i) {
      for (size_t j=0; j<base.size(); ++j) {
         cout << base[i] << " x " << base[j] << " score="
            << nuc.lookup(base[i], base[j]) << " ";
      }
      cout << endl;
   }

   nuc.use("NUC.4.4.N");
   cout << "min score: " << nuc.getMinScore() << endl;
   cout << "max score: " << nuc.getMaxScore() << endl
      << " gap open: " << nuc.getGapOpen() << " gap extend: " << nuc.getGapExtend()
      << endl;
   //nuc.show();
   //string base("ACGTN");
   for (size_t i=0; i<base.size(); ++i) {
      for (size_t j=0; j<base.size(); ++j) {
         cout << base[i] << " x " << base[j] << " score="
            << nuc.lookup(base[i], base[j]) << " ";
      }
      cout << endl;
   }

   cout << "test constructor\n";
   NucleicScoreMethod nucN("NUC.4.4.N");
   cout << "min score: " << nucN.getMinScore() << endl;
   cout << "max score: " << nucN.getMaxScore() << endl
      << " gap open: " << nucN.getGapOpen() << " gap extend: " << nucN.getGapExtend()
      << endl;
   nucN.show();
   for (size_t i=0; i<base.size(); ++i) {
      for (size_t j=0; j<base.size(); ++j) {
         cout << base[i] << " x " << base[j] << " score="
            << nucN.lookup(base[i], base[j]) << " ";
      }
      cout << endl;
   }
   cout << "end of test nucleic\n";
}

int main(int argc, char *argv[]) {
   MatrixScoreMethod::setDefaultPath("../matrix");
   testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
