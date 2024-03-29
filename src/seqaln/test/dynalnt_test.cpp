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
      DynalnTest() : smbase(), ssm(10, -11, -40, -1),
         psm(), nsm(), nsmN("NUC.4.4.N"),
         alnsmbase(smbase),
         alnssm(ssm),
         alnpsm(psm),
         alnnsm(nsm)
        { }

      virtual void SetUp() { }
      virtual void TearDown() { }

      void testlocal(const DNA& sq1, const DNA& sq2) {
         alnssm.setSeq(sq1, sq2);
         alnssm.runlocal();
         alnssm.printAlign(cout, 80);
      }

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

TEST_F(DynalnTest, printaln) {
   SimpleScoreMethod sm(10, -11, -40, -1);
   DNA seq1("sq1", "GTCCCCTCACAGTTCTGAGGTCTGAAATCAAGGTGTCTGCAGAGCTGGTTCCTGTTGAGGGTTGCAAAGGAGAATCTCTTCCAGCCTCTCTCCTAGCTTTGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTATACATGTAGATATATGTGTGCATGTGCATATATGTGCTCACATGTGTACATGC");
   DNA seq2("sq2", "GTCCCCTCACAGTTCTGAGGTCTGAAATCAAGGTGTCTGCAGAGCTGGTTCCTGTTGAGGGTTGCAAAGCAGAATCTCTTCCAGCCTCTCTCCTAGCTTTGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTAAACATGTAGATATATGTGT");
   Dynaln<SimpleScoreMethod> aligner(sm);
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   aligner.printAlign(cout);
}

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
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   ASSERT_TRUE(aligner.getIdentity() > 0.5);
   string nonseq = aligner.getNCorrectedBottomSequence();
   cout << "N-corrected\n" << nonseq << endl
      << "original:\n" << nseq << endl;
   ASSERT_TRUE(nonseq.size() == nseq.length());
   seq1.setSequence("AGAATCACTGGCAGTAGTGGACAGCACCTGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGATAAGCATGTCACAAATACATTCTCTGATCTAATCCTCACCACACACACAATGGGCACTATTCCTTCTACAGATGAGGAAACTGA");
   nseq.setSequence("GGCCAGAACCTGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTNTGTGATAAGCATGTCACAAATACATTCTCTGATCTAATCCTCACCACACACACAATGGGCACTATTCCTTCTACAGATGAGGAAA");
   aligner.setseq(seq1, nseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   ASSERT_TRUE(aligner.getIdentity() > 0.5);
   nonseq = aligner.getNCorrectedBottomSequence();
   cout << "N-corrected\n" << nonseq << endl
      << "original:\n" << nseq << endl;
   ASSERT_TRUE(nonseq.size() == nseq.length());

   cout << "Another test make sure no N left in the sequence\n";
   // test N still left or not
   seq1.setSequence("aatcatctaatggaatggaatggaataatccatggactcgaatgcaatcatcatcgaatggaatcgaatggaatcatcgaatggactcgaatggaataatcattgaacggaatcgaatggaatcatcatcggatggaaatgaatggaatcatcatcgaatggaatcgaatagaattatgga"); // lower case
   nseq.setSequence("AATGGAATGGAAAAGAATGGACTGTNANTCAANAATCATCGAANGGAATNGANTGGAATCNTCGAATGGACTCGAATGGAATAATCATTGAACGGAATCGAATGGAATCATCATCGGATGGAAATGAATGGAATCATCATCGAATGGAATN");
   aligner.setseq(seq1, nseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   nonseq = aligner.getNCorrectedBottomSequence();
   cout << "N-corrected\n" << nonseq << endl
      << "original:\n" << nseq << endl;
   ASSERT_TRUE(nonseq.find('N') == string::npos);
}

TEST_F(DynalnTest, fixstagger) {
   DNA seq1("topseq", "GCTTTAGGAATCGGATACTTAGCCAAATTAAAAGGCAAAAGATAAACCAGGAAATACACTTTTAATAAATAAGCCAAAAGATAACTGTTTTTGATATGTAAAGAAGTCAATAGGGGTGGGGCACAGTGGCTCATGCCTATAATTCCAGCACTTTGGGAGGCTGAGACGGGCGGATCTCTTG");
   DNA qseq("queryseq", "ATTCCTTGGAAGCAGGACAGGAATTGAGTCATATCCACTCATGGCACCAGGTGCCATAATAGGCGCTCAATAAATGTTTGCTGTGGCAGGGCGCGGTGGCTCATGCCTATAATTCCAGCAC");
// 81        91        101       111       118       128       138       148
// +         +         +         +         +         +         +         +         
// ATAACTGTTTTTGATATGTAAAGAAGTCAATA---GGGGTGGGGCACAGTGGCTCATGCCTATAATTCCAGCAC
// |||| |   |||                       | ||  |||| | ||||||||||||||||||||||||||
// ATAAAT--GTTT--------------------GCTGTGGCAGGGCGCGGTGGCTCATGCCTATAATTCCAGCAC
// +         +         +         +         +         +         +         +  
// 70        78                            88        98        108       118

   SimpleScoreMethod ssm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethod> aligner(ssm);
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   ASSERT_TRUE(aligner.fixStagger());
   cout << "after fix stagger\n";
   aligner.printAlign(cout, 80);
   vector<pair<char,int>> cigar=aligner.getCigar1();
   cout << "Cigar:\n";
   for (auto& p : cigar) {
      cout << p.second << p.first;
   }
   cout << endl;
   cout << aligner.getMDString1() << endl;
   // test zag configuration
//102                           124       134       144       154       164
//+         +         +         +         +         +         +         +
//GAGGGATGG--GGA-GGA-----CATGACCCCCCGAGCCACCTTCCCTGCCGGGCCTTTCCAGCCGTCCCAGAGCCAGTC
//|| ||| ||  ||| |||          ||||||||||||||||||| |||||||||| ||||||||||| ||||| |||
//GACGGACGGACGGACGGACGGGG-----CCCCCCGAGCCACCTTCCCCGCCGGGCCTTCCCAGCCGTCCCGGAGCCGGTC
//+         +         +         +         +         +         +         +
//69        79        89        94        104       114       124       134

   seq1.setSequence("CCCGGGATTCTGCGAGTGCTACTGCTGAGGGACTGTAACACTCGGGGTGTGGCCCAGCTCACCCCCCTCCAGAGGGATGGGGAGGACATGACCCCCCGAGCCACCTTCCCTGCCGGGCCTTTCCAGCCGTCCCAGAGCCAGTCACGGCGCAGCACCACAGTGGAAATGCATCTGGCGGTGG");
   qseq.setSequence("CGCCCCGACCCGCGCGCCCTCCCGCGGGAGGACGCGGGGCCGGGGGGCGGAGACGGGGGAGGAGGACGGACGGACGGACGGACGGACGGGGCCCCCCGAGCCACCTTCCCCGCCGGGCCTTCCCAGCCGTCCCGGAGCCGGTCGCGGCGCA");
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   //cerr << "MD before fixStagger: " << aligner.getMDString1() << endl;
   ASSERT_TRUE(aligner.fixStagger());
   cout << "after fix stagger\n";
   aligner.printAlign(cout, 80);
   cigar=aligner.getCigar1();
   cout << "Cigar:\n";
   for (auto& p : cigar) {
      cout << p.second << p.first;
   }
   cout << endl;
   cout << aligner.getMDString1() << endl;
}

TEST_F(DynalnTest, fixStaggerGap) {
   DNA seq1("topseq", "GCTTTAGGAATCGGATACTTAGCCAAATTAAAAGGCAAAAGATAAACCAGGAAATACACTTTTAATAAATAAGCCAAAAGATAACTGTTTTTGATATGTAAAGAAGTCAATAGGGGTGGGGCACAGTGGCTCATGCCTATAATTCCAGCACTTTGGGAGGCTGAGACGGGCGGATCTCTTG");
   DNA qseq("queryseq", "ATTCCTTGGAAGCAGGACAGGAATTGAGTCATATCCACTCATGGCACCAGGTGCCATAATAGGCGCTCAATAAATGTTTGCTGTGGCAGGGCGCGGTGGCTCATGCCTATAATTCCAGCAC");
// 81        91        101       111       118       128       138       148
// +         +         +         +         +         +         +         +         
// ATAACTGTTTTTGATATGTAAAGAAGTCAATA---GGGGTGGGGCACAGTGGCTCATGCCTATAATTCCAGCAC
// |||| |   |||                       | ||  |||| | ||||||||||||||||||||||||||
// ATAAAT--GTTT--------------------GCTGTGGCAGGGCGCGGTGGCTCATGCCTATAATTCCAGCAC
// +         +         +         +         +         +         +         +  
// 70        78                            88        98        108       118

   SimpleScoreMethod ssm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethod> aligner(ssm);
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   ASSERT_TRUE(aligner.fixStaggerGap());
   cout << "after fix stagger\n";
   aligner.printAlign(cout, 80);
   vector<pair<char,int>> cigar=aligner.getCigar1();
   cout << "Cigar:\n";
   for (auto& p : cigar) {
      cout << p.second << p.first;
   }
   cout << endl;
   cout << aligner.getMDString1() << endl;
   // test zag configuration
//102                           124       134       144       154       164
//+         +         +         +         +         +         +         +
//GAGGGATGG--GGA-GGA-----CATGACCCCCCGAGCCACCTTCCCTGCCGGGCCTTTCCAGCCGTCCCAGAGCCAGTC
//|| ||| ||  ||| |||          ||||||||||||||||||| |||||||||| ||||||||||| ||||| |||
//GACGGACGGACGGACGGACGGGG-----CCCCCCGAGCCACCTTCCCCGCCGGGCCTTCCCAGCCGTCCCGGAGCCGGTC
//+         +         +         +         +         +         +         +
//69        79        89        94        104       114       124       134

   seq1.setSequence("CCCGGGATTCTGCGAGTGCTACTGCTGAGGGACTGTAACACTCGGGGTGTGGCCCAGCTCACCCCCCTCCAGAGGGATGGGGAGGACATGACCCCCCGAGCCACCTTCCCTGCCGGGCCTTTCCAGCCGTCCCAGAGCCAGTCACGGCGCAGCACCACAGTGGAAATGCATCTGGCGGTGG");
   qseq.setSequence("CGCCCCGACCCGCGCGCCCTCCCGCGGGAGGACGCGGGGCCGGGGGGCGGAGACGGGGGAGGAGGACGGACGGACGGACGGACGGACGGGGCCCCCCGAGCCACCTTCCCCGCCGGGCCTTCCCAGCCGTCCCGGAGCCGGTCGCGGCGCA");
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.fixStaggerGap());
   cout << "after fix stagger gap\n";
   aligner.printAlign(cout, 80);
   cigar=aligner.getCigar1();
   cout << "Cigar:\n";
   for (auto& p : cigar) {
      cout << p.second << p.first;
   }
   cout << endl;
   cout << aligner.getMDString1() << endl;
   /*
    * TCTTC--AGGTGATATTCCAGATACCATG
    * |||||   |||||||||||||||||||||
    * TCTTCGT-GGTGATATTCCAGATACCATG
    */
   seq1.setSequence("CAAGCCATTTCCAGCCATTTGATTTTCTTCAGGTGATATTCCAGATACCATGAAACAGTGACAAGCCATCCCCATGAAGTCCTGTCTAAATTCTTGACCCAT");
   qseq.setSequence("GATCTTCGTGGTGATATTCCAGATACCATGAAACAGTGACAAGCCATCCCCATGAAGTCCTGTCTAAATTCTTGACCCAT");
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.fixStaggerGap());
   cout << "after fix stagger gap\n";
   aligner.printAlign(cout, 80);
   /* 
44        54        64        74        82        92        102       112      
+         +         +         +         +         +         +         +         
TGGAAAGGAATGGAATAGAATGGAATCTTCCTGAGA--GGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGG
|| |||||||||||| |||||  ||||    ||| |   | |||||||||| ||||||||||||||||||||||||||||
TGAAAAGGAATGGAAGAGAATAAAATCCAG-TGAAATT-GGATGGAATGGATTGGAATGGAATGGAATGGAATGGAATGG
+         +         +         +         +         +         +         +         
5         15        25        31        43        53        63        73       
After fix generates another stagger

44        54        64        74        83        93        103       113
+         +         +         +         +         +         +         +
TGGAAAGGAATGGAATAGAATGGAATCTTCCTGAGA-GGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGA
|| |||||||||||| |||||  ||||      |    | |||||||||| |||||||||||||||||||||||||||||
TGAAAAGGAATGGAAGAGAATAAAATCCAGTGAAATT-GGATGGAATGGATTGGAATGGAATGGAATGGAATGGAATGGA
+         +         +         +         +         +         +         +
5         15        25        35        44        54        64        74

*/
   seq1.setSequence("TGGAAAGGAATGGAATAGAATGGAATCTTCCTGAGAGGAATGGAATGGAATGGAATGGAATGGAATGGAATGGAATGG");
   qseq.setSequence("TGAAAAGGAATGGAAGAGAATAAAATCCAGTGAAATTGGATGGAATGGATTGGAATGGAATGGAATGGAATGGAATGG");
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.fixStaggerGap());
   cout << "after fix stagger gap\n";
   aligner.printAlign(cout, 80);
   /**
    * Fixing algorithm crashed
177       187       197       207       217       227       237       247
+         +         +         +         +         +         +         +
GAGTGGAATGGAATGGAACGAAATGGAATAGAATGGAACGGAATGCAATGGAATGGATGGAAATTGAATGGAAAGGAATA
|||||||||||||||||| |||||||||                     ||||||||      | |||||     |||
GAGTGGAATGGAATGGAATGAAATGGAA--------------------AGGAATGGAA-----TAGAATG-----GAA--
+         +         +         +         +         +         +         +
1         11        21        31        41        31        61        71

257       267       275       285       295       305
+         +         +         +         +         + 
AAAACAAGTGAAAT--TGGATGGAATGGATTGGAATGGAATGGAATGGAAT
             |  | || |||||||| |||||||||||||||||||||
----------TCTTCCTAGA-GGAATGGAATGGAATGGAATGGAATGGAAT
+         +         +         +         +         + 
81        49        101       68        78        88

*/
   seq1.setSequence("GAGTGGAATGGAATGGAACGAAATGGAATAGAATGGAACGGAATGCAATGGAATGGATGGAAATTGAATGGAAAGGAATAAAAACAAGTGAAATTGGATGGAATGGATTGGAATGGAATGGAATGGAAT");
   qseq.setSequence("GAGTGGAATGGAATGGAATGAAATGGAAAGGAATGGAATAGAATGGAATCTTCCTAGAGGAATGGAATGGAATGGAATGGAATGGAAT");
   aligner.setseq(seq1, qseq);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   aligner.fix1M();
   cout << "after fix1M()\n";
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.fixStaggerGap());
   cout << "after fix stagger gap\n";
   aligner.printAlign(cout, 80);

}


TEST_F(DynalnTest, fix1M) {
   DNA seq1("refseq", "AAAATGTCTAAGCAAGTCAAAGAGCATTTATGAATAATAGGCCTGTGAGAAAACTTTTATGAATGATCAGG");
   DNA seq2("queryseq", "AAAATTAATAACAAATTAAAATATCATTTATCAAAAAAAAAAAAAATAAAAAACTTTTATGAATGATCAGG");
   SimpleScoreMethodN ssm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethodN> aligner(ssm);
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   bool fixed=aligner.fix1M();
   if (fixed) {
      cout << "after fixing 1M problem\n";
      aligner.printAlign(cout, 80);
   }
   else {
      ASSERT_TRUE(false);
   }
   aligner.fixStagger();
   cout << endl << "after fixing stagger\n";
   aligner.printAlign(cout, 80);
   aligner.fixStaggerGap();
   cout << endl << "after fixing stagger gap\n";
   aligner.printAlign(cout, 80);
   /*
CTTATTTATT--T-ATAGCTCCTATCATTATTATTTCCTCCCTCTTATAACTTTAGGT
||| || |||    | |||||||||||||||| |||||||||||||||||||||||||
CTTTTTAATTGA-AACAGCTCCTATCATTATTTTTTCCTCCCTCTTATAACTTTAGGT
*/
   seq1.setSequence("TGACAGTTTTGCCTGTTTTAATTAGCTGTTTCAAATAACAAGCTCCTGGCACATTTATTTAAATTTCATGACATGTTTTCTTATTTATTTATAGCTCCTATCATTATTATTTCCTCCCTCTTATAACTTTAGGTTTGTTGTTTTTGCTGTTCTTTTTCATATGTGAAATGGAATGGTTTTTTTGATAATGTTGACATTTTTAAA");
   seq2.setSequence("AAGAAATTTACATAAATGTGCCAGCAGCTTTTTAATTGAAACAGCTCCTATCATTATTTTTTCCTCCCTCTTATAACTTTAGGT");
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   aligner.printAlign(cout, 80);
   aligner.fix1M();
   cout << endl << "after fixing 1M\n";
   aligner.printAlign(cout, 80);
   /*
1         11        21        27        36        45        55        65        75
+         +         +         +         +         +         +         +         +
TTTTGTTTTTCCCTCTCATTCTCCTT----ACAGCTC-TACCCACT-CCTTCTTTCTTCAACAGATATTTACTGAGTATTTT
||||||||||   ||| |||||  ||    | ||    ||   |     || |||  |  | ||||||| | ||||||||||
TTTTGTTTTTAGATCTGATTCTAATTTAATATAGAAAATAT--A--A--TTATTTAGTAGAAAGATATTGAATGAGTATTTT
+         +         +         +         +         +         +         +         +
1         11        21        31        41        45        55        65        75
*/
   seq1.setSequence("TTTTGTTTTTCCCTCTCATTCTCCTTACAGCTCTACCCACTCCTTCTTTCTTCAACAGATATTTACTGAGTATTTT");
   seq2.setSequence("TTTTGTTTTTAGATCTGATTCTAATTTAATATAGAAAATATAATTATTTAGTAGAAAGATATTGAATGAGTATTTT");
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   aligner.printAlign(cout, 90);
   //ASSERT_TRUE(aligner.fixStaggerGap());
   ASSERT_TRUE(aligner.fix1M());
   cout << "after fix 1M\n";
   aligner.printAlign(cout, 90);
   /*
1         11        21        27        36        46        56        66        76
+         +         +         +         +         +         +         +         +
TTTTGTTTTTCCCTCTCATTCTCCTT----ACAGCTC-TACCCACTCCTTCTTTCTTCAACAGATATTTACTGAGTATTTT
||||||||||   ||| |||||  ||    | ||    ||        || |||  |  | ||||||| | ||||||||||
TTTTGTTTTTAGATCTGATTCTAATTTAATATAGAAAATATA-----ATTATTTAGTAGAAAGATATTGAATGAGTATTTT
+         +         +         +         +         +         +         +         +
1         11        21        31        41        46        56        66        76
*/
   ASSERT_TRUE(aligner.fixStaggerGap());
   cout << "after fix stagger gap\n";
   aligner.printAlign(cout, 90);
}

TEST_F(DynalnTest, trimLeft) {
   DNA seq1("refseq", "TGTTAGCTTCATATTAAAGCTGTTTTGTTTCATGGTCAAAAGATAGATGCCAATGGTGGCAATGGATGCCTAGCAGGACAGAGACTCACTACACACACAAACACACACACACACACACACACACACACACACACACACACACAAAGTAGAT");
   DNA seq2("query", "TTTTTTATATATATTTATTATCTTTTTTTTCTTTTTTTCAAAATAGATATCACCCTTGTGTCTCTATTTCTCGCTCTAGAGAGAGAGACACTAAACACACACACACACACACACACACACACACACACACACACACACACACACCCTTTTC");
   /*
ref x S41733017  10-147/148 | 10-150/151
Score=576 gap length: 5 2 num gaps: 3 1 idencnt=106 simcnt=0 alnlen=143 identity=0.74126 similarity=0.74126
11        21        31        39        49        59        69        77
+         +         +         +         +         +         +         +
ATATTAAAGCTGTTTTGTTTCATGGT--CAAAAGATAGATGCCAATGGTGGCAATGGATGCCTAG--CAGGACAGAGACT
||||| |   | |||| |||| |  |  ||||| | | ||  |  | ||| |  |  ||  || |  |  || |||||
ATATTTATTATCTTTTTTTTCTTTTTTTCAAAATAGATATCACCCTTGTGTC--TCTATTTCTCGCTCTAGAGAGAGAGA
+         +         +         +         +         +         +         +
11        21        31        41        51        61        69        79

87        97        107       117       127       136       146
+         +         +         +         +         +         +         +
CACTACACACACAAACACACACACACACACACACACACACACACACAC-CACACACCCTTTTC
||||| ||||||| |||||||||||||||||||||||||||||||||| ||||||||||||||
CACTAAACACACACACACACACACACACACACACACACACACACACACACACACACCCTTTTC
+         +         +         +         +         +         +         +
89        99        109       119       129       139       149
*/
   SimpleScoreMethod ssm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethod> aligner(ssm);
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   cout << "Before trimming\n";
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   aligner.trimLeft(0.85);
   cout << "after left trimming. gaplen2=" << aligner.getGaplen2() << "\n";
   aligner.printAlign(cout, 80);
   //cout << "after print  gaplen2=" << aligner.getGaplen2() << "\n";
   ASSERT_TRUE(aligner.getAlnlen() == 56);
   // another more difficult case
   seq1.setSequence("ACCACAGTGGAATTAAATTAGAATTACTCACCAAAAGTTAACTAGGGGAAACACAAACACACACACACACACACACACACACACACACAGAGCGACACAGACAAGATAGCGCGAGCGCACGGGTCAGCTCGCCTCTCTCGGCTGGAGCTCG");
   seq2.setSequence("TATACTTTGGGTATAACTCTTTTTCTACACTCCGCTCTCTAGATATCTCACAGAGTGTTTATATATTATATTATCTCAAAAAATTATATCTCGGTTACACACACACACACACACACACACACACACACACACACAGAGAGCCACATACACA");
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   cout << "Before trimming\n";
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   aligner.trimLeft(0.85);
   cout << "after left trimming. gaplen2=" << aligner.getGaplen2() << "\n";
   aligner.printAlign(cout, 80);
   // TTTGAATACTAGAGTGTTATCGTTTTTTTTTTATTTTTTTTTAATTTTTTTTTTTTT
   // ||| |||| || | | || |   ||||||||||||||||||||||||||||||||||
   // TTTTAATATTA-A-TATTTTTT-TTTTTTTTTATTTTTTTTTAATTTTTTTTTTTTT
   seq1.setSequence("TTTGAATACTAGAGTGTTATCGTTTTTTTTTTATTTTTTTTTAATTTTTTTTTTTTT");
   seq2.setSequence("TTAGTTTTTGGTAGGGATCAATTTTAATATTAATATTTTTTTTTTTTTTTATTTTTTTTTAATTTTTTTTTTTTTATCAATTTTTATTTTTTTAATAATATATTAATTACTATAAAATAATGTTAATTAAATACTAAATCTTTTTATAAAA");
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   cout << "Before trimming\n";
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   aligner.trimLeft(0.85);
   cout << "after left trimming. gaplen2=" << aligner.getGaplen2() << "\n";
   aligner.printAlign(cout, 80);
   //ref: AGCGGGAGGGGGCGGGCAGGGACACTTACACGCTCGCCAGGGGGTCCGGGCAGGCCAGTGGGTCCGGGTTCACCTCGAGCTCCTCGCTTTCCTTCCTCTCCAGCAGGGCCAGGGGGTCCTTGAACACCAACAGGGCC
   //query: CCGATCTACCTCCCCAGGTTTGCCTGCTTCACCTGGAGGACCAGCAGGTCCAGGGAGACCCTGGAAGCCGGGGGAGCCAGCAGGGCCTTGTTCACCTCTCTCGCCAGCGGGACCAGCAGGGCCAGGGGGTCCCTGAACACCAACAGGGCCAGGCG
   seq1.setSequence("AGCGGGAGGGGGCGGGCAGGGACACTTACACGCTCGCCAGGGGGTCCGGGCAGGCCAGTGGGTCCGGGTTCACCTCGAGCTCCTCGCTTTCCTTCCTCTCCAGCAGGGCCAGGGGGTCCTTGAACACCAACAGGGCC");
   seq2.setSequence("CCGATCTACCTCCCCAGGTTTGCCTGCTTCACCTGGAGGACCAGCAGGTCCAGGGAGACCCTGGAAGCCGGGGGAGCCAGCAGGGCCTTGTTCACCTCTCTCGCCAGCGGGACCAGCAGGGCCAGGGGGTCCCTGAACACCAACAGGGCCAGGCG");
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   cout << "Before fixStagger\n";
   aligner.printAlign(cout, 80);
/* before fix stagger
refseq x query  36-136/137 | 49-149/155
Score=434 gap length: 18 18 num gaps: 3 4 idencnt=74 simcnt=0 alnlen=119 identity=0.62185 similarity=0.62185
37        45        44        57        67        77        87        97
+         +         +         +         +         +         +         +
CCAGGGGG--TCCGGGCAG--------GCCAGTGGGTCCGGGTTCACCTCGAGCTCCTCGCTTTCCTTCCTCTC------
|||||| |   || || ||        |||||  || ||  |||||||||    | |||||            |
CCAGGGAGAC-CCTGGAAGCCGGGGGAGCCAGCAGGGCCTTGTTCACCTC----T-CTCGC------------CAGCGGG
+         +         +         +         +         +         +         +
50        52        69        79        89                  104

          109       119       129
+         +         +         +         +         +         +         +
--CAGCAGGGCCAGGGGGTCCTTGAACACCAACAGGGCC
  ||||||||||||||||||| |||||||||||||||||
ACCAGCAGGGCCAGGGGGTCCCTGAACACCAACAGGGCC
+         +         +         +         +         +         +         +
112       122       132       142      
*/
   aligner.fixStagger();
   cout << "Before trimming, after fixStagger\n";
/*
refseq x query  36-136/137 | 49-149/155
Score=434 gap length: 9 9 num gaps: 2 3 idencnt=76 simcnt=0 alnlen=110 identity=0.69091 similarity=0.69091
37        46        44        58        68        78        88        98
+         +         +         +         +         +         +         +
CCAGGGGG-TCCGGGCAG--------GCCAGTGGGTCCGGGTTCACCTCGAGCTCCTCGCTTTCCTTCCTCTCCAGCAGG
|||||| |  || || ||        |||||  || ||  |||||||||    | |||||    |  |    ||||||||
CCAGGGAGACCCTGGAAGCCGGGGGAGCCAGCAGGGCCTTGTTCACCTC----T-CTCGC----CAGCGGGACCAGCAGG
+         +         +         +         +         +         +         +         
50        60        70        80        90                  104       111

108       118       128       129
+         +         +         +         +         +         +         +
GCCAGGGGGTCCTTGAACACCAACAGGGCC
|||||||||||| |||||||||||||||||
GCCAGGGGGTCCCTGAACACCAACAGGGCC
+         +         +         +         +         +         +         +
121       131       141       142
*/
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.getNumgaps1() == 2 && aligner.getNumgaps2() == 3);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   // testing trim left
   aligner.trimLeft(0.85);
   cout << "after left trimming. gaplen2=" << aligner.getGaplen2() << "\n";
/*
refseq x query  94-136/137 | 107-149/155
Score=434 gap length: 0 0 num gaps: 0 0 idencnt=38 simcnt=0 alnlen=43 identity=0.88372 similarity=0.88372
95        105       115       125       135       7
+         +         +         +         +         +         +         +
CCTCTCCAGCAGGGCCAGGGGGTCCTTGAACACCAACAGGGCC
|    |||||||||||||||||||| |||||||||||||||||
CGGGACCAGCAGGGCCAGGGGGTCCCTGAACACCAACAGGGCC
+         +         +         +         +         +         +         +
108       118       128       138       148
*/
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.getAlnlen() == 43);
   /*
refseq x query  82-228/299 | 14-158/159
Score=659 gap length: 2 4 num gaps: 1 2 idencnt=110 simcnt=0 alnlen=149 identity=0.73826 similarity=0.73826
83        93        103       113       123       133       143       153
+         +         +         +         +         +         +         +
CCCGAGAGGTGGAGGTTGCAGTGAGCCAAGATCACACCACTGCACTCTCGTCTAGGTGATAGAGCGAGACTCTGTCTC--
|| | | || ||||  ||| |||||||  | |  | |||| || | |  | |  || |  |   | | || |     |
CCGGGGGGGGGGAGCCTGCTGTGAGCCGCGGTGGCGCCACGGCGCCCCAGCC--GGGGCCA--CCCAAACCCCCAAACAA
+         +         +         +         +         +         +         +
15        25        35        45        55        65        73        81

161       171       181       191       201       211       221
+         +         +         +         +         +         +         +
AAAAAAAAAAAAAAAAAAAAAAAAAGTAACTCCAAGGCTGGGCGTGGTGGCTCATGCCTGTAATCCCAG
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| | ||
AAAAAAAAAAAAAAAAAAAAAAAAAGTAACTCCAAGGCTGGGCGTGGTGGCTCATGCCTGTAATACAAG
+         +         +         +         +         +         +         +
91        101       111       121       131       141       151
*/
   seq1.setSequence("AAATACAAAAATTAGCTGGGGGTGGTGGCGGGTGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGGAGCAGAATCGCTTGAACCCGAGAGGTGGAGGTTGCAGTGAGCCAAGATCACACCACTGCACTCTCGTCTAGGTGATAGAGCGAGACTCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAGTAACTCCAAGGCTGGGCGTGGTGGCTCATGCCTGTAATCCCAGCACATTAGGAGGCCGAAGTAGGCAGATTGCTTGAGGCCAGGAGTTCGAGACCAGCCTGGCCAACATAGTG");
   seq2.setSequence("GCTCTTCCTATCCCCCGGGGGGGGGGAGCCTGCTGTGAGCCGCGGTGGCGCCACGGCGCCCCAGCCGGGGCCACCCAAACCCCCAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAGTAACTCCAAGGCTGGGCGTGGTGGCTCATGCCTGTAATACAAG");
   testlocal(seq1, seq2);
   alnssm.trimLeft(0.91);
   cout << "after trim left at 0.91\n";
   alnssm.printAlign(cout, 80);
}

TEST_F(DynalnTest, trimRight) {
   DNA seq1("refseq", "TCATGACTGACAAGGGGAGGGGAAGTAACACTTTTTAAGTTCACCCTGGGAAAGTAACCACTATCTACACACACACCAGTACATACACACACACACACACACACACACACACACCAGTACACACACACACACACACACACACACACACACACCAGTACA");
   DNA seq2("query", "TCATGACTGACAAGGGGAGGGGAAGTAACACTTTTTAAGTTCACCCTGGGAAAGTAACCACTATCTACACACACACCAGTACATACACACACACACACACACAATAACAACAGCAAAACACAAAACACAACAAACAAGGACAGATCTAGTG");
/*
1         11        21        31        41        51        61        71
+         +         +         +         +         +         +         +
TCATGACTGACAAGGGGAGGGGAAGTAACACTTTTTAAGTTCACCCTGGGAAAGTAACCACTATCTACACACACACCAGT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
TCATGACTGACAAGGGGAGGGGAAGTAACACTTTTTAAGTTCACCCTGGGAAAGTAACCACTATCTACACACACACCAGT
+         +         +         +         +         +         +         +
1         11        21        31        41        51        61        71

81        91        101       111       121       130       140       150
+         +         +         +         +         +         +         +
ACATACACACACACACACACACACACACACACACCAGTACACACACA-CACACACACACACACACACACACACCAGTACA
|||||||||||||||||||||||            | || ||| ||| || | ||||| | ||||| ||| || || |||
ACATACACACACACACACACACA------------A-TA-ACA-ACAGCA-AAACACAAA-ACACA-ACAAACAAGGACA
+         +         +         +         +         +         +         +
81        91        101                 107                           133
*/
   SimpleScoreMethod ssm(10, -11, -40, -1);
   Dynaln<SimpleScoreMethod> aligner(ssm);
   aligner.setseq(seq1, seq2);
   aligner.runlocal();
   cout << "Before trimming\n";
   aligner.printAlign(cout, 80);
   cout << "nogap identity=" << aligner.getNogapIdentity() << endl;
   aligner.trimRight(0.85);
   cout << "after right trimming\n";
   aligner.printAlign(cout, 80);
   ASSERT_TRUE(aligner.getAlnlen() == 103);
   // another test
   /*
ref x S29252  70-240/299 | 0-147/159
Score=674 gap length: 5 28 num gaps: 2 8 idencnt=122 simcnt=0 alnlen=176 identity=0.69318 similarity=0.69318
71        81        91        101       111       121       131       141
+         +         +         +         +         +         +         +
CGCCCGGGGCCGCTGCAGATGGCGGGGCGCGTTGGAGCGCCGGGCCGGCCCCGGGGCTGGAGGGAGGCCCGCGAGACCCC
|||||||||||||||||||||||||||||||||||||||||||||||| ||||||    |  || |  || | |  ||||
CGCCCGGGGCCGCTGCAGATGGCGGGGCGCGTTGGAGCGCCGGGCCGGGCCCGGG---CGCCGGCG--CCACCACGCCCC
+         +         +         +         +         +         +         +
1         11        21        31        41        51        58        66

151       169       166       176       186       196       206       216
+         +         +         +         +         +         +         +
GGGCCGTC----CGCCC-CGCCGCCGCGCTCCGGCCCGCGGGGGCAGCTTGCGAGCCCCGACGCCCCGGGCCCAGGGCCG
||| || |    ||    ||||||||      ||| ||| | |||     |  ||| |||     ||||||||     ||
GGGGCGGCGGGGCGGGGGCGCCGCCG------GGCGCGCCGAGGC-----GG-AGCTCCG-----CCGGGCCC-----CG
+         +         +         +         +         +         +         +
76        86        96        119       110       115                 129

226       236
+         +         +         +         +         +         +         +
CGCTCCGAAGCGCCGC
||| |||  |||||||
CGC-CCGGTGCGCCGC
+         +         +         +         +         +         +         +
134       143
*/
   seq1.setSequence("CCTGCGGAGCTACGCCTCGGAGCGCCCGTCGGCGGCCCCGACCCGCAGTCCCCGGGCCTGGATGAGCCTGCGCCCGGGGCCGCTGCAGATGGCGGGGCGCGTTGGAGCGCCGGGCCGGCCCCGGGGCTGGAGGGAGGCCCGCGAGACCCCGGGCCGTCCGCCCCGCCGCCGCGCTCCGGCCCGCGGGGGCAGCTTGCGAGCCCCGACGCCCCGGGCCCAGGGCCGCGCTCCGAAGCGCCGCTTCCAGAACTCGACCCGTTGTTCTCCTGGACTGAGGAGCCCGAGGAGTGTGGCCCCGC");
   seq2.setSequence("CGCCCGGGGCCGCTGCAGATGGCGGGGCGCGTTGGAGCGCCGGGCCGGGCCCGGGCGCCGGCGCCACCACGCCCCGGGGCGGCGGGGCGGGGGCGCCGCCGGGCGCGCCGAGGCGGAGCTCCGCCGGGCCCCGCGCCCGGTGCGCCGCCCCCCGCCCGG");
   testlocal(seq1,seq2);
   alnssm.trimRight(0.92);
   cout << "After trim right at 0.92\n";
   alnssm.printAlign(cout, 80);
   /*
    * >ref
GGTGAACCTCAGCACGTCATTGACGTCGGTGGAGGACACGGCGCTCCACCCTGGTGGCCGAGGGTGGGCTGCCTCCCCTC
CTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTCATTCAGCAAATATTTATTTATTTATTTTG
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTTATTGTTTATTTTTTATTTCTCCAAACCTCTTTTTTTTTTTATTGATCATT
CTTGGGTGTTTCTCACAGAGGGGGATTTGGCAGGGTCATAGGACAATAGTGGAG

>S11232167
GCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTCATTCAGCAAATATTTATTTAT
TTATTGTGTGTGTGTGTGTGTGTGTGTGTGTACAGACAGAATTTTTTTTCACTGTTTTTTTATATGTATCGA
ref x S11232167  70-236/294 | 0-149/152
Score=1075 gap length: 2 19 num gaps: 2 5 idencnt=135 simcnt=0 alnlen=169 identity=0.79882 similarity=0.79882
71        81        91        101       111       121       131       141
+         +         +         +         +         +         +         +
GCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTCATTCAGCAAATATTTATTT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||  ||||||||||||||||||||||
GCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT--CTCATTCAGCAAATATTTATTT
+         +         +         +         +         +         +         +
1         11        21        31        41        51        59        69

151       161       171       181                 200       209       219
+         +         +         +         +         +         +         +
ATTTATTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTT-ATTGTTTATTTT-TTATTTCTCCAAACCTCTTTTTTTTT
|||||||    ||||||||||||||||||||||||||    |  |         || ||| |      | || |||||||
ATTTATT----GTGTGTGTGTGTGTGTGTGTGTGTGTACAGACAG-------AATTTTTT-TT-----CACTGTTTTTTT
+         +         +         +         +         +         +         +
79                  95        105       115                           132

229
+         +         +         +         +         +         +         +
TTATTGATC
 |||  |||
ATATGTATC
+         +         +         +         +         +         +         +
142

Trim is bad: 
ref x S11232167  70-190/294 | 0-115/152
Score=1075 gap length: 1 6 num gaps: 1 2 idencnt=112 simcnt=0 alnlen=122 identity=0.91803 similarity=0.91803
71        81        91        101       111       121       131       141
+         +         +         +         +         +         +         +
GCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTCATTCAGCAAATATTTATTT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||  ||||||||||||||||||||||
GCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT--CTCATTCAGCAAATATTTATTT
+         +         +         +         +         +         +         +
1         11        21        31        41        51        59        69

151       161       171       181
+         +         +         +         +         +         +         +
ATTTATTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTT-A
|||||||    ||||||||||||||||||||||||||    |
ATTTATT----GTGTGTGTGTGTGTGTGTGTGTGTGTACAGA
+         +         +         +         +         +         +         +
79                  95        105       115
*/
   seq1.setSequence("GGTGAACCTCAGCACGTCATTGACGTCGGTGGAGGACACGGCGCTCCACCCTGGTGGCCGAGGGTGGGCTGCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTCATTCAGCAAATATTTATTTATTTATTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTATTGTTTATTTTTTATTTCTCCAAACCTCTTTTTTTTTTTATTGATCATTCTTGGGTGTTTCTCACAGAGGGGGATTTGGCAGGGTCATAGGACAATAGTGGAG");
   seq2.setSequence("GCCTCCCCTCCTGTGATCAGCAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTCATTCAGCAAATATTTATTTATTTATTGTGTGTGTGTGTGTGTGTGTGTGTGTACAGACAGAATTTTTTTTCACTGTTTTTTTATATGTATCGA");
   testlocal(seq1,seq2);
   alnssm.trimRight(0.8);
   cout << "After trim right at 0.8\n";
   alnssm.printAlign(cout, 80);
   cout << string(20, '-') << " trimRight() " << string(20, '-') << endl;
   /*
Score=379 gap length: 0 0 num gaps: 0 0 idencnt=72 simcnt=0 alnlen=103 identity=0.69903 similarity=0.69903
21        31        41        51        61        71        81        91       
+         +         +         +         +         +         +         +         
TGTGAGCCACTGTGCCAGGCCTGGCCTTTTCTTTTACTCATTGGACTTCGTGGTGACCCAATCCCGTTTCTTCCCCTGTT
|||||||||||||||||||  |||  |||| |||   | ||||  ||| ||| | |   |||   |||| ||   || | 
TGTGAGCCACTGTGCCAGGAATGGAATTTTATTTGTATAATTGCTCTTTGTGCTCAAAAAATAATGTTTATTAAACTCTG
+         +         +         +         +         +         +         +         
1         11        21        31        41        51        61        71       

101       111       121       1
+         +         +         +         +         +         +         +         
AGAGGGCACATTTAGAATTACTT
| |||  ||||| | ||||||||
ACAGGCTACATTGACAATTACTT
+         +         +         +         +         +         +         +         
81        91        101       1
     test trimRight(0.91)
*/
   Dynaln<SimpleScoreMethod>::setTrimWidth(26);
   seq1.setSequence("TGTGAGCCACTGTGCCAGGCCTGGCCTTTTCTTTTACTCATTGGACTTCGTGGTGACCCAATCCCGTTTCTTCCCCTGTTAGAGGGCACATTTAGAATTACTT");
   seq2.setSequence("TGTGAGCCACTGTGCCAGGAATGGAATTTTATTTGTATAATTGCTCTTTGTGCTCAAAAAATAATGTTTATTAAACTCTGACAGGCTACATTGACAATTACTT");
   testlocal(seq1,seq2);
   alnssm.printAlign(cout,80);
   cout << "above is before trimRight()\n";
   alnssm.trimRight(0.89);
   cout << "After trim right at 0.89\n";
   alnssm.printAlign(cout, 80);
}
