#include "scorematrix.h"
#include "bioseq.h"
#include "dynalnt.h"
#include "strformat.h"
#include <gtest/gtest.h>

using namespace std;
using namespace orpara;

/**
 * For testing alignments. 
 * This is a fixture named pattern: ClassNameTest
 */
class DynalnTest : public testing::Test {
   protected:
      DynalnTest() : smbase(), ssm(20, -19),
         psm(), nsm(), nsmN("NUC.4.4.N"),
         alnsmbase(smbase),
         alnssm(ssm),
         alnpsm(psm),
         alnnsm(nsm)
        { }

      virtual void SetUp() { }
      virtual void TearDown() { }

      ScoreMethod smbase;
      SimpleScoreMethod ssm;
      ProteinScoreMethod psm;
      NucleicScoreMethod nsm;
      NucleicScoreMethod nsmN;

      Dynaln<ScoreMethod> alnsmbase;
      Dynaln<SimpleScoreMethod> alnssm;
      Dynaln<ProteinScoreMethod> alnpsm;
      Dynaln<NucleicScoreMethod> alnnsm;
};

TEST_F(DynalnTest, repeataln) {
   string sA1 = "ACCCGTGTGGCCAGTCGTACGTGTACACTGACGTACGATCAGTC";
   string sB1 = "ACGTGTCCAGTCCGTGCTGACGAATTGTGTGTTTTGGTTTC";
   string sA2 = "ACGGGTGTGCCCAGACGT";
   string sB2 = "TTTTTTGGGGTTTGGC";

   SimpleScoreMethod sm(5, -4, -8, -8);
   int maxLen = max(max(sA1.size(), sA2.size()), max(sB1.size(), sB2.size()));

   cout << "seqA : " << sA1 << " -- " << sA2 << endl;
   cout << "seqB : " << sB1 << " -- " << sB2 << endl;

   Dynaln<SimpleScoreMethod> aln(sm);
   //double normalScore[maxLen + 1], score12, score11, score22;
   double score12, score11, score22;

   //for (int l = 1; l <= maxLen; l++) {
   //for (int l = 3; l <= maxLen; l++) {
   for (int l = 3; l <= 5; l++) {
      int x, y;
      x=max<int>(0, int(sA1.size())-l);
      y=max<int>(0, int(sB1.size()) - l);
      if (x > int(sA1.size()-1) || y > (int)sB1.size()) {
         throw runtime_error("x, y passed end of input string, cannot take substring beond end");
      }
      cout << string(2 * maxLen + 20, '=') << endl;
      cout << "\nl=" << l << " x y: " << x << " " << y << endl;

      DNA parSeqA1("seqA1", sA1.substr(x, l));
      DNA parSeqB1("seqB1", sB1.substr(y, l));
      DNA parSeqA2("seqA2", sA2.substr(0, l));
      DNA parSeqB2("seqB2", sB2.substr(0, l));

      aln.setSeq(parSeqA1, parSeqB1);
      aln.runlocal();
      score11 = aln.getScore();
      cout << parSeqA1 << parSeqB1 << endl;
      aln.printAlign(cout);
      //aln.debug_showmatrix(cout);
      cout << string(30, '&') << endl;

      //cout << "parSeqA1: " << parSeqA1.toString() << endl;
      //cout << "parSeqB2: " << parSeqB2.toString() << endl;
      aln.setSeq2(parSeqB2);
      aln.runlocal();
      score12 = aln.getScore();
      cout << parSeqA1 << parSeqB2 << endl;
      aln.printAlign(cout);
      //aln.debug_showmatrix(cout);
      cout << string(30, '&') << endl;

      //cout << "parSeqA2: " << parSeqA2.toString() << endl;
      //cout << "parSeqB2: " << parSeqB2.toString() << endl;
      aln.setSeq1(parSeqA2);
      aln.runlocal();
      score22 = aln.getScore();
      cout << parSeqA2 << parSeqB2 << endl;
      aln.printAlign(cout);
      //aln.debug_showmatrix(cout);
      cout << string(30, '&') << endl;
      cout << "score 1x1=" << score11 << " 1x2=" << score12 << " 2x2=" << score22 << endl;
   }
}

TEST_F(DynalnTest, simplescoremethodshortseq) {
   SimpleScoreMethod sm(5, -4, -8, -8);
   Dynaln<SimpleScoreMethod> aln;
   double tmp1=0, tmp2;

   DNA seqA2("seqA2", "ACG");
   DNA seqB2("seqB2", "TTT");
   aln.setseq(seqA2, seqB2);
   aln.runlocal();
   tmp1=aln.getScore();
   aln.printAlign(cout);
   cout << "transition to length 4\n";

   DNA seqA1("seqA1", "AGTC");
   DNA seqB1("seqB1", "TTTC");
   aln.setseq(seqA1, seqB1);
   aln.runlocal();
   tmp2=aln.getScore();
   aln.printAlign(cout);
   cout << tmp1 << " " << tmp2 << endl;
   ASSERT_GT(aln.getScore(), 7);

   DNA dna1("sq1", "ACGGGTG");
   DNA dna2("sq2", "TTTTTTG");
   aln.setseq(dna1, dna2);
   aln.runlocal();
   cout << "Aligning two sequences:\n" << dna1 << dna2 << endl;
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 7);
   DNA dna3("seq3", "CGATCAGTC");
   DNA dna4("seq4", "TTTGGTTTC");
   aln.setseq(dna3, dna4);
   aln.runlocal();
   aln.printAlign(cout);
}

TEST_F(DynalnTest, nucmatrix) {
   DNA seq1("CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGACGAAACACCGG");
   DNA seq2("CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGACAACACCGG");
   cout << "test score matrix base\n";
   alnsmbase.setSeq(seq1, seq2);
   alnsmbase.runlocal();
   alnsmbase.printAlign(cout);
   ASSERT_GT(alnsmbase.getScore(), 0);

   cout << "test score matrix smple\n";
   alnssm.setSeq(seq1, seq2);
   alnssm.runlocal();
   alnssm.printAlign(cout);
   ASSERT_GT(alnssm.getScore(), 0);

   cout << "test score matrix default NUC.4.4\n";
   alnnsm.setSeq(seq1, seq2);
   alnnsm.runlocal();
   alnnsm.printAlign(cout);
   ASSERT_GT(alnnsm.getScore(), 0);
   string base("ACGTN");
   cout << "NUC.4.4 base similarity\n";
   for (size_t i=0; i<base.size(); ++i) {
      for (size_t j=0; j<base.size(); ++j) {
         cout << base[i] << "x" << base[j] << "="
            << nsm.lookup(base[i], base[j]) << "  ";
      }
      cout << endl;
   }

   cout << "test score matrix NUC.4.4.N\n";
   NucleicScoreMethod nucmatrixN("NUC.4.4.N");
   cout << "default path: " << NucleicScoreMethod::getDefaultPath() << endl;
   cout << "gap open: " << nucmatrixN.getGapOpen()
      << " gap extend: " << nucmatrixN.getGapExtend()
      << endl;
   cout << "NUC.4.4.N base similarity\n";
   for (size_t i=0; i<base.size(); ++i) {
      for (size_t j=0; j<base.size(); ++j) {
         cout << base[i] << "x" << base[j] << "="
            << nucmatrixN.lookup(base[i], base[j]) << "  ";
      }
      cout << endl;
   }
   Dynaln<NucleicScoreMethod> naln(nucmatrixN);
   //naln.setMatrix(nucmatrixN);
   naln.setSeq(seq1, seq2);
   naln.runlocal();
   naln.printAlign(cout);
   ASSERT_GT(naln.getScore(), 0);
}

TEST_F(DynalnTest, prtaln1) {
   Protein prt1("prt1", "IPYEPEPELLTTPASAPSQ");
   Protein prt2("prt2", "IQYDPDPELLSTPASSPSN");
   /* take two input files for alignment input */
   Dynaln<ProteinScoreMethod> aln2prt(prt1, prt2);
   aln2prt.runlocal();
   aln2prt.printAlign(cout);
   ASSERT_GT(aln2prt.getScore(), 0);

   Protein p1("Protein1", "IPYEPEPELLTTTPASAPSQSTSRKRKLDEPDSRPSTNRIVELDSVPLDLSRAAADDLRRRVAAQSRMLDASERWLAAVASNAPGTTDDPSARHLADSCQSRVDALHKLSHASPLRLHFASNGARVAQSLAKRQ");
   Protein p2("Protein2", "TPYEPEPELLAERASSKTPTEAFDDATNRKRKHEALATAIPEASRFEEIDSAPLDNVRAVLDDLRRDVDAHARFADASERWLAALATNTPGTAEDPESFADAGESDAHLTAFNASVAAERAAEARLRLLLDAAP");
   Dynaln<ProteinScoreMethod> aln(p1,p2);
   aln.runglobal();
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 0);
}

TEST_F(DynalnTest, proteinAlign) {
   //Protein prt1("P28329", "MAAKTPSSEESGLPKLPVPPLQQTLATYLQCMRHLVSEEQFRKSQAIVQQFGAPGGLGETLQQKLLERQEKTANWVSEYWLNDMYLNNRLALPVNSSPAV");
   //Protein prt2("NP_066265", "MAAKTPSSEESGLPKLPVPPLQQTLATYLQCMRHLVSEEQFRKSQAIVQQFGAPGGLGETLQQKLLERQEKTANWVSEYWLNDMYLNNRLALPVNSSPAV");
   /* expected alignment
   prt1 1   MTGIEHVVKRYYIDGSVYALPHDQGEHSRLTKLHKLFLRILDHKILHAPV     50
            |:||||..||.||||.|||||||..|.|||...|.||||..|:|:|||||
   prt2 138 MSGIEHKAKRDYIDGVVYALPHDTEESSRLLMQHNLFLRAFDNKVLHAPV    187
   */
   Protein prt1("prt1", "MTGIEHVVKRYYIDGSVYALPHDQGEHSRLTKLHKLFLRILDHKILHAPV");
   Protein prt2("prt2", "MSGIEHKAKRDYIDGVVYALPHDTEESSRLLMQHNLFLRAFDNKVLHAPV");
   cout << prt1 << prt2 << endl;
   ProteinScoreMethod defaultsm;
   defaultsm.setGapOpen(-32);
   Dynaln<ProteinScoreMethod> aln2prt(defaultsm);
   aln2prt.setseq(prt1, prt2);
   aln2prt.runlocal();
   //aln2prt.runglobal();
   aln2prt.printAlign(cout);
   string top = aln2prt.getTopAln();
   string middle = aln2prt.getMiddleAln();
   string bottom = aln2prt.getBottomAln();
   EXPECT_STREQ(top.c_str(), "MTGIEHVVKRYYIDGSVYALPHDQGEHSRLTKLHKLFLRILDHKILHAPV");
   EXPECT_STREQ(middle.c_str(), "|:|||| :|| |||| |||||||  | |||   |:|||| :|:|:|||||");
   EXPECT_STREQ(bottom.c_str(), "MSGIEHKAKRDYIDGVVYALPHDTEESSRLLMQHNLFLRAFDNKVLHAPV");
}

TEST_F(DynalnTest, proteinMixedCase) {
   cout << "testing mixed case sequence for protein alignment\n";
   Protein prt1("P28329", "MGLRTAKKRGLGGGGKWKREEGGGTRGRREVRPACFLQSGGRGDPGDVGGPAGNPGCSPHPRAATRPPPLPAHTPAHTPEWCGAASAEAAEPRRAGPHLCIPAPGLTKTPILEKVPRKMAAKTPSSEESGLPKLPVPPLQQTLATYLQCMRHLVSEEQFRKSQAIVQQFGAPGGLGETLQQKLLERQEKTANWVSEYWLNDMYLNNRLALPVNSSPAVIFARQHFPGTDDQLRFAASLISGVLSYKALLDSHSIPTDCAKGQLSGQPLCMKQYYGLFSSYRLPGHTQDTLVAQNSSIMPEPEHVIVACCNQFFVLDVVINFRRLSEGDLFTQLRKIVKMASNEDERLPPIGLLTSDGRSEWAEARTVLVKDSTNRDSLDMIERCICLVCLDAPGGVELSDTHRALQLLHGGGYSKNGANRWYDKSLQFVVGRDGTCGVVCEHSPFDGIVLVQCTEHLLKHVTQSSRKLIRADSVSELPAPRRLRWKCSPEIQGHLASSAEKLQRIVKNLDFIVYKFDNYGKTFIKKQKCSPDAFIQVALQLAFYRLHRRLVPTYESASIRRFQEGRVDNIRSATPEALAFVRAVTDHKAAVPASEKLLLLKDAIRAQTAYTVMAITGMAIDNHLLALRELARAMCKELPEMFMDETYLMSNRFVLSTSQVPTTTEMFCCYGPVVPNGYGACYNPQPETILFCISSFHSCKETSSSKFAKAVEESLIDMRDLCSLLPPTESKPLATKEKATRPSQGHQP");
   string pepstr = "maaktpsseesglpklpvpplqqtlatylqcmrhlvseeqfrksqaivqqfgapgglgetlqqkllerqektanwvseywlndmylnnrlalpvnsspavifarqhfpgtddqlrfaaslisgvlsykalldshsiptdcakgqlsgqplcmkqyyglfssyrlpghtqdtlvaqnssimpepehvivaccnqffvldvvinfrrlsegdlftqlrkivkmasnederlppiglltsdgrsewaeartvlvkdstnrdsldmierciclvcldapggvelsdthralqllhgggysknganrwydkslqfvvgrdgtcgvvcehspfdgivlvqctehllkhmtqssrkliradsvselpaprrlrwkcspeiqghlassaeklqrivknldfivykfdnygktfikkqkcspdafiqvalqlafyrlhrrlvptyesasirrfqegrvdnirsatpealafvravtdhkaavpasekllllkdairaqtaytvmaitgmaidnhllalrelaramckelpemfmdetylmsnrfvlstsqvptttemfccygpvvpngygacynpqpetilfcissfhscketssskfakaveeslidmrdlcsllppteskplatkekatrpsqghqp";

   //Protein prt2("NP_066265", str2upper(pepstr));
   Protein prt2("NP_066265", pepstr);
   cout << prt1 << prt2 << endl;

   ProteinScoreMethod defaultsm;
   //Dynaln<ProteinScoreMethod> aln2prt(prt1, prt2);
   defaultsm.setGapOpen(-12);
   Dynaln<ProteinScoreMethod> aln2prt(defaultsm);
   aln2prt.setseq(prt1, prt2);
   aln2prt.runlocal();
   //aln2prt.runglobal();
   aln2prt.printAlign(cout);
   ASSERT_GT(aln2prt.getScore(), 0);
}

TEST_F(DynalnTest, badNucGlobal) {
   cerr << "Testing bad input DNA sequece with global algorithm\n";
   // not sure what to expect this bad input
   DNA d1("top", "ACCAAAAAACAACCCCCCCCCCCC");
   DNA d2("bottom", "GGGTTTTTTTGGGGGGGGGGGGTT");
   SimpleScoreMethod sm(10, -9, -11, -3);
   Dynaln<SimpleScoreMethod> alnr(d1, d2, sm);
   alnr.runglobal();
   alnr.printAlign(cout);
   ASSERT_LT(alnr.getScore(), 0);

   cout << "chang to a good set of sequence\n";
   DNA g1("top", "GCTCTCTTTTCAGAGCAACAAACTC");
   DNA g2("top", "GCTCTCTTTCAGAGCAAACAAACTC");
   alnr.setSeq(g1, g2);
   alnr.showseq();
   alnr.runglobal();
   alnr.printAlign(cout);
   alnr.showseq();
   ASSERT_GT(alnr.getScore(), 0);
}

TEST_F(DynalnTest, testBadNucLocal) {
   cerr << "Testing bad input DNA sequece with local algorithm\n";
   DNA d1("top", "ACCAAAAAACAACCCCCCCCCCCC");
   DNA d2("bottom", "GGGTTTTTTTGGGGGGGGGGGGTT");
   SimpleScoreMethod sm(10, -9, -11, -3);
   Dynaln<SimpleScoreMethod> alnr(d1, d2, sm);
   alnr.runlocal();
   alnr.printAlign(cout);
   ASSERT_EQ(alnr.getScore(), 0);

   cout << "chang to a good set of sequence\n";
   DNA g1("top", "GCTCTCTTTTCAGAGCAACAAACTC");
   DNA g2("top", "GCTCTCTTTCAGAGCAAACAAACTC");
   alnr.setSeq(g1, g2);
   alnr.showseq();
   alnr.runlocal();
   alnr.printAlign(cout);
   alnr.showseq();
   ASSERT_GT(alnr.getScore(), 0);
}

TEST_F(DynalnTest, linearSpaceProtein) {
   cout << "\nLinear space testing global protein alignment\n";
   Protein p1("protein1", "IPYEPEPELLTTTPASAPSQSTSR");
   Protein p2("protein2", "TPYEPEPELLAERASSKTPTEA");
   LSDynaln<ProteinScoreMethod> aln1(p1,p2);
   aln1.global();
   aln1.buildResult();
   aln1.printAlign(cout);

   // not testing all the others
   cout << endl << "Testing Local Linear protein alignment\n";
   aln1.local();
   aln1.buildResult();
   aln1.printAlign(cout);
   ASSERT_GT(aln1.getScore(), 0);
}


TEST_F(DynalnTest, nucleicAcidAlign) {
   cerr << "Testing nucleic acid alignment ...\n";
   /*
   Matrix scm;
   scm.setDefaultPath();
   scm.setMatrix("NUC.4.4");
   */
   SimpleScoreMethod spcm(5, -6);
   spcm.show();

   DNA seq1("CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGACGAAACACCGG");
   DNA seq2("CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGACAACACCGG");

   //Dynaln<NucleicScoreMethod> aln;
   Dynaln<SimpleScoreMethod> aln;
   aln.setMatrix(spcm); // this will fail because Simple cannot be convert to
   // NucleicScoreMethod
   aln.setGapInsert(-9);
   aln.setGapExtend(-2);
   spcm.show();

   aln.setSeq1(seq1);
   aln.setSeq2(seq2);
   aln.runlocal();
   aln.printAlign(cout);

   // template
   cout << "sequence assigned from string object\n";
   seq1 = "TTTTTGAATTCTCGACCTCGAGACAAATGGCAGTATTCATCCACAAGATCGGAAGAGC";
   seq2 = "CCTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGACGAAACACCGGCCTAACTGTGCATGGCAACTCGAGTTGCCATGACAACCAGTTAGGGTTTTGAATTTCTCGACCTCGAGACAAATGGCAGTATTCATCCACAGATCGGAAGAGC";
   seq1.show();
   seq2.show();
   aln.setSeq1(seq1);
   aln.setSeq2(seq2);
   aln.runlocal();
   aln.printAlign(cout);

   seq1 = "GGCGTATATGACAACAATATCTC";
   seq2 = "ATATGACAACCAATATCTC";
   aln.setSeq1(seq1);
   aln.setSeq2(seq2);
   aln.runglobal();
   aln.printAlign(cout);
   cout << "no terminal gap identity: " << aln.getNoTerminalGapIdentity() << endl;
   ASSERT_GT(aln.getScore(), 0);
}

TEST_F(DynalnTest, nucMatrix) {
   cout << "Testing nucleic acid matrix 4x4 ...\n";
   DNA dna1("site", "GACGGATCGGGAGATCTCCCGATCCCCTATGGTGCACTCTCAGTACAATCTG");
   DNA dna2("read", "GACGGATCGGAGATCTCCCGATCCCCTATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCAT");
   //dna1.getcode();
   //dna2.getcode();
   dna1.show();
   dna2.show();
   //MatrixScoreMethod::setDefaultPath("/home/zhouke/src/proj/seqaln/matrix");
   NucleicScoreMethod scm("NUC.4.4");
   //scm.show();
   cout << "G=G " << scm.lookup('G', 'G')
      << " A=A " << scm.lookup('A', 'A')
      << " C=C " << scm.lookup('C', 'C')
      << " T=T " << scm.lookup('T', 'T')
      << " A=G " << scm.lookup('A', 'G') << endl;
   Dynaln<NucleicScoreMethod> aln;
   aln.setMatrix(scm);
   aln.setGapParameter(-6, -4);
   //cout << "showing all parameters:\n";
   //aln.showParameters();
   //cout << "parameters shown\n";
   cout << "Local Matrix 4.4. alignment\n";
   aln.setSeq1(dna1);
   aln.setSeq2(dna2);
   aln.runlocal();
   aln.printAlign(cout);
   cout << "NUC.4.4 matrix test done\n";
   ASSERT_GT(aln.getScore(), 0);

   // now use identity matrix
   cout << "Using identity matrix ...\n";
   SimpleScoreMethod idscm(6, -5);
   Dynaln<SimpleScoreMethod> alni(dna1, dna2);
   alni.setMatrix(idscm);
   alni.setGapInsert(-4);
   alni.setGapExtend(-4);
   alni.runlocal();
   alni.printAlign(cout);
   cout << "identity matrix test success\n";
   ASSERT_GT(alni.getScore(), 0);
}

TEST_F(DynalnTest, longGap) {
   /*
    * When the back trace pointer was not computed properly,
    * the alignment may be wrong.
    * Now the pointer calculation is a more
    * complicated and using two piece of information.
    * It will give the correct alignment.
    */
   cout << "\n\ntesting long gap\n";
   DNA d1("shorseq", "gccaattccgcccctctccccactcccctctcgcca");
   DNA d2("longseq", "gccaattccgcccctctccccaccacccctctccctcctaacgttactgtcccctctcgcca");
   SimpleScoreMethod scm(5, -4);
   Dynaln<SimpleScoreMethod> aln;
   aln.setMatrix(scm);
   aln.setGapParameter(-30, -1);
   aln.setSeq1(d1);
   aln.setSeq2(d2);
   aln.runlocal();
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 0);
}

TEST_F(DynalnTest, doubleGaps) {
   /*
    * For some iotorrent sequences, gaps are frequent. We need to set the
    * proper gap score relative to the match score. We only need the simple
    * scoring scheme for matching sequence to itself. This is only modeling
    * the machine error, no need for biological concern.
    */
   string expectedAln=
         string("AGGGAGTCTGGCGGGGCGACGGCATC ATGCAAACCACCTGCCCATGTGGAGCACAGATCA\n")
      +  string("|||||||||||||||||||||||||| |||||| ||||||||||||||||||||| |||||\n")
      +  string("AGGGAGTCTGGCGGGGCGACGGCATCCATGCAA CCACCTGCCCATGTGGAGCAC-GATCA\n");

   cout << "\n\ntesting double gap\n"
      << "expected alignment:\n" << expectedAln << endl;
   DNA d1("refseq", "AGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAGATCA");
   DNA d2("read1", "AGGGAGTCTGGCGGGGCGACGGCATCCATGCAACCACCTGCCCATGTGGAGCACGATCA");
   SimpleScoreMethod scm(10, -9, -30, -3);
   Dynaln<SimpleScoreMethod> aln;
   aln.setMatrix(scm);
   aln.setSeq(d1, d2);
   aln.runlocal();
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 0);

   expectedAln=
        string("CCGGACATGTGAAAAACGGTTCC-ATGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAA\n")
      + "||||||||||||||| ||||||| |||||||||||||| |||||||||||||||||||||||||||||||\n"
      + "CCGGACATGTGAAAA-CGGTTCCTATGAGGATCGTGGG-CCTAGGACCTGTAGTAACACGTGGCATGGAA\n";

   d1.setSequence("CCGGACATGTGAAAAACGGTTCCATGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAA");
   d2.setSequence("CCGGACATGTGAAAACGGTTCCTATGAGGATCGTGGGCCTAGGACCTGTAGTAACACGTGGCATGGAA");

   cout << "expected alignment:\n" << expectedAln << endl;
   aln.setSeq(d1, d2);
   aln.runlocal();
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 0);

   d1.setSequence("GTACGTGGAGGTTACGCGGGTGGGGGATTTCCATT");
   d2.setSequence("GTACGTGGAGGTTACGCGGCAGGGGGGATTTCCACT");
   aln.setSeq(d1, d2);
   aln.setGapInsert(-30);
   aln.runlocal();
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 0);

   d1.setSequence("ACGGTTCCATGAGGATCGTTGGGCCTAAAACCTGCAGCAACACGTGGCATGGAACATTCCCCAT");
   d2.setSequence("ACGGTTCCATGAGGATCGTTGGGCCAAAGACCTGCAGCAACACGTGGCATGGAACATTCCCAT");
   aln.setSeq(d1, d2);
   aln.setGapInsert(-29);
   aln.runlocal();
   aln.printAlign(cout);
   ASSERT_GT(aln.getScore(), 0);
}

TEST_F(DynalnTest, nucmatrixN) {
   DNA seq1("i5adapterRC", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT");
   DNA seq2("read_cause_segfault", "CGCTCAAGGCTCACCCCTCCTGCCCTGTGTCCCTACAGACACTAACAGCACATCTGGAGACCCGGTGGAGAAGAAGGACGAAACAGTTAGATCGGAAGAGCGTCGAGTAGGGAAAGAGTGTTTGGACTGGTGTAGATCTCGGGGGTCGCCG");
   cout << "test score matrix with N\n";
   Dynaln<NucleicScoreMethod> aligner(nsmN);

   aligner.setSeq(seq1, seq2);
   aligner.runlocal();
   aligner.printAlign(cout);
   ASSERT_GT(aligner.getScore(), 0);
}

TEST_F(DynalnTest, simplematrixLongN) {
   DNA seq1("goodseq", "GTACAAGGCACTGTGCTATATCTGGAACTACAAATATTAGTTAAATAGTCCAAGATTCAGTGGGTTCAGAGACTCTAAGTATTTGCCATTCCCTTTAGATGGAAACGTTCTTTGGTCCACCAAGCCCTATCGTTCTGCTCTTTGGAAATTCTCCCTCACTTTCTTCATAATGCCGTCCTTGACTAACTCTCCTCTAGCAATAC");
   DNA nseq("seqwithN", "CAAATATTAGNTAANNNNTCCAAGATTNANTGNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTAGATGGAAACGTTCTTTGGTCCACCANGCCCTATCGNNNTGCTCNTTGGNNATTCTNCCTCACNTTCTTCNTAATGC");
   SimpleScoreMethodN ssm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethodN> aligner(ssm);
   aligner.setseq(seq1, nseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.getIdentity() > 0.5);
}

