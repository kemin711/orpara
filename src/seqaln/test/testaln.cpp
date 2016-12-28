#include "scorematrix.h"
#include "bioseq.h"
#include "dynalnt.h"
#include "strformat.h"

void testNucleicAcidAlign();
void testNucMatrix();
void testLongGap();
bool testLinearSpaceProtein();
void testBadNuc();
void testBadNucLocal();
void testProteinMixedCase();
void testProteinAlign();
void testDoubleGaps();
void testSpeed();

/** 
 * a program to tes the alignemnt classes.
 */
int main(int argc, char *argv[]) {
   Protein prt1, prt2;
   /* take two input files for alignment input */
   if (argc > 1) {
      prt1.read(argv[1]);
      prt2.read(argv[2]);
      Dynaln<ProteinScoreMethod> aln2prt(prt1, prt2);
      aln2prt.runlocal();
      aln2prt.printAlign(cout);
      return 0;
   }

   Protein p1("Protein1", "IPYEPEPELLTTTPASAPSQSTSRKRKLDEPDSRPSTNRIVELDSVPLDLSRAAADDLRRRVAAQSRMLDASERWLAAVASNAPGTTDDPSARHLADSCQSRVDALHKLSHASPLRLHFASNGARVAQSLAKRQ");
   Protein p2("Protein2", "TPYEPEPELLAERASSKTPTEAFDDATNRKRKHEALATAIPEASRFEEIDSAPLDNVRAVLDDLRRDVDAHARFADASERWLAALATNTPGTAEDPESFADAGESDAHLTAFNASVAAERAAEARLRLLLDAAP");
   Dynaln<ProteinScoreMethod> aln(p1,p2);
   aln.runglobal();
   aln.printAlign(cout);


   /*
   cerr << "local testing\n";
   p1.assign("IPYEPEPELLTTTPASAPSQSTSR");
   p2.assign("TPYEPEPELLAERASSKTPTEA");
   aln.setSeq1(p1);
   aln.runlocal();
   //aln.buildResult();
   aln.printAlign(cout);
   */

   // when both functions are used, I got segmentation fault!
   testNucleicAcidAlign();
   testNucMatrix();
   testLongGap();
   testBadNucLocal();
   testBadNuc();
   // the following is not working
   testLinearSpaceProtein();
   testProteinMixedCase();
   testProteinAlign();
   testDoubleGaps();

   return 0;
}
void testProteinAlign() {
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
}

void testProteinMixedCase() {
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
}


void testBadNuc() {
   cerr << "Testing bad input DNA sequece with global algorithm\n";
   DNA d1("top", "ACCAAAAAACAACCCCCCCCCCCC");
   DNA d2("bottom", "GGGTTTTTTTGGGGGGGGGGGGTT");
   SimpleScoreMethod sm(10, -9, -11, -3);
   Dynaln<SimpleScoreMethod> alnr(d1, d2, sm);

   alnr.runglobal();
   alnr.printAlign(cout);

   cout << "chang to a good set of sequence\n";
   DNA g1("top", "GCTCTCTTTTCAGAGCAACAAACTC");
   DNA g2("top", "GCTCTCTTTCAGAGCAAACAAACTC");
   alnr.setSeq(g1, g2);
   alnr.showseq();
   alnr.runglobal();
   alnr.printAlign(cout);
   alnr.showseq();
}

void testBadNucLocal() {
   cerr << "Testing bad input DNA sequece with local algorithm\n";
   DNA d1("top", "ACCAAAAAACAACCCCCCCCCCCC");
   DNA d2("bottom", "GGGTTTTTTTGGGGGGGGGGGGTT");
   SimpleScoreMethod sm(10, -9, -11, -3);
   Dynaln<SimpleScoreMethod> alnr(d1, d2, sm);

   alnr.runlocal();
   alnr.printAlign(cout);

   cout << "chang to a good set of sequence\n";
   DNA g1("top", "GCTCTCTTTTCAGAGCAACAAACTC");
   DNA g2("top", "GCTCTCTTTCAGAGCAAACAAACTC");
   alnr.setSeq(g1, g2);
   alnr.showseq();
   alnr.runlocal();
   alnr.printAlign(cout);
   alnr.showseq();
}

bool testLinearSpaceProtein() {
   cout << "\nLinear space testing global protein alignment\n";
   Protein p1("protein1", "IPYEPEPELLTTTPASAPSQSTSR");
   Protein p2("protein2", "TPYEPEPELLAERASSKTPTEA");
   LSDynaln<ProteinScoreMethod> aln1(p1,p2);
   aln1.global();
   aln1.buildResult();
   aln1.printAlign(cout);

   //* not testing all the others
   cout << endl << "Testing Local Linear protein alignment\n";
   aln1.local();
   aln1.buildResult();
   aln1.printAlign(cout);
   return true;
}


void testNucleicAcidAlign() {
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

}

void testNucMatrix() {
   cout << "Testing nucleic acid matrix 4x4 ...\n";
   DNA dna1("site", "GACGGATCGGGAGATCTCCCGATCCCCTATGGTGCACTCTCAGTACAATCTG");
   DNA dna2("read", "GACGGATCGGAGATCTCCCGATCCCCTATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCAT");
   //dna1.getcode();
   //dna2.getcode();
   dna1.show();
   dna2.show();
   MatrixScoreMethod::setDefaultPath("/home/zhouke/src/proj/seqaln/matrix");
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
}

void testLongGap() {
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
}

void testDoubleGaps() {
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

   d1.setSequence("GTACGTGGAGGTTACGCGGGTGGGGGATTTCCATT");
   d2.setSequence("GTACGTGGAGGTTACGCGGCAGGGGGGATTTCCACT");
   aln.setSeq(d1, d2);
   aln.setGapInsert(-30);
   aln.runlocal();
   aln.printAlign(cout);

   d1.setSequence("ACGGTTCCATGAGGATCGTTGGGCCTAAAACCTGCAGCAACACGTGGCATGGAACATTCCCCAT");
   d2.setSequence("ACGGTTCCATGAGGATCGTTGGGCCAAAGACCTGCAGCAACACGTGGCATGGAACATTCCCAT");
   aln.setSeq(d1, d2);
   aln.setGapInsert(-29);
   aln.runlocal();
   aln.printAlign(cout);


}

/** give a pair of sequences test the speed */
void testSpeed() {
}
