#include "scorematrix.h"

#include <vector>
#include "strformat.h"
#include <iomanip>

///// class ScoreMethod  the base class //////////////
namespace orpara {
void ScoreMethod::show(ostream &ous) const {
   ous << "ScoreMethod class with gapOpen=" << gapo << " gapExtend=" << gape << endl;
}

ScoreMethod& ScoreMethod::operator=(const ScoreMethod &scm) {
   if (this != &scm) {
      gapo = scm.gapo;
      gape = scm.gape;
   }
   return *this;
}

ScoreMethod& ScoreMethod::operator=(ScoreMethod &&scm) {
   if (this != &scm) {
      gapo = scm.gapo;
      gape = scm.gape;
      scm.gapo=0;
      scm.gape=0;
   }
   return *this;
}


/////// class SimpleScoreMethod ////////////////////////
SimpleScoreMethod& SimpleScoreMethod::operator=(const SimpleScoreMethod &mt) {
   if (this != &mt) {
      matchs = mt.matchs;
      mismatchs = mt.mismatchs;
      setGapParameters(mt.gapo, mt.gape);
   }
   return *this;
}


////////// MatrixScoreMethod   the the scoring matrix class ///////////

///// static method section /////////////
//

//string MatrixScoreMethod::default_path="/home/zhouke/src/proj/seqaln/matrix";
//string MatrixScoreMethod::default_path="/usr/local/share/orpara/matrix";
// MATRIXDIR will be replaced by the autotool automake
// AM_CPPFLAGS = -D MATRIXDIR="\"$(pkgdatadir)/matrix\""
string MatrixScoreMethod::default_path=MATRIXDIR;
const char* MatrixScoreMethod::proteinMatrices[]  = {
   // put all new matrices here
   "dayhoff", "gonnet", "identity", "match", 
   "blosumn", "blosumn.50", 
   "blosum30", "blosum30.50", "blosum35", "blosum35.50", "blosum40", "blosum40.50", 
   "blosum45", "blosum45.50", "blosum50", "blosum50.50", "blosum55", "blosum55.50", 
   "blosum60", "blosum60.50", "blosum62", "blosum62.50", "blosum65", "blosum65.50", 
   "blosum70", "blosum70.50", "blosum75", "blosum75.50", "blosum80", "blosum80.50", 
   "blosum85", "blosum85.50", "blosum90", "blosum90.50", "blosum100", "blosum100.50", 
   "pam10", "pam20", "pam30", "pam40", "pam40.cdi", "pam50", "pam60", "pam70", 
   "pam80", "pam80.cdi", "pam90" "pam100", "pam110", "pam120", "pam120.cdi", 
   "pam130", "pam140", "pam150", "pam160", "pam160.cdi", "pam170", "pam180", 
   "pam190", "pam200", "pam200.cdi", "pam210", "pam220", "pam230", "pam240", 
   "pam250", "pam250.cdi", "pam260", "pam270", "pam280", "pam290", "pam300", 
   "pam310", "pam320", "pam330", "pam340", "pam350", "pam360", "pam370", "pam380", 
   "pam390", "pam400", "pam410", "pam420", "pam430", "pam440", "pam450", "pam460", 
   "pam470", "pam480", "pam490", "pam500"
};

// when you add elements to this array, the number 
// in isNucleicMatrix also needs to be incremented!
const char* MatrixScoreMethod::nucleicMatrices[] = {
   "DNA", "DNAelem", "NUC.4.4", "NUC.4.4.N"
};

bool MatrixScoreMethod::isProteinMatrix(const string &matrixName) {
   string lowername = getLower(matrixName);
   if (lowername.substr(0,6) == "blosum"
         || lowername.substr(0,3) == "pam") {
      return true;
   }
   else {
      for (int i = 0; i<4; ++i) {
         if (proteinMatrices[i] == lowername) return true;
      }
   }
   return false;
}

// may need to add more logic
bool MatrixScoreMethod::isNucleicMatrix(const string &matrixName) {
   for (int i = 0; i<4; ++i) {
      if (nucleicMatrices[i] == matrixName) return true;
   }
   return false;
}

void MatrixScoreMethod::setDefaultPath() {
   char *matrixDir = getenv("SCORING_MATRIX_PATH");
   if (matrixDir != NULL) {
      default_path = matrixDir;
   }
   else {
      //cerr << "SCORING_MATRIX_PATH environment variable is not set, and I will look for the scoring matrix in the the source directory\n";
      //string pathToMatrix(getenv("HOME"));
      //pathToMatrix += "/src/proj/seqaln/matrix";
      //default_path=pathToMatrix;
      default_path="/usr/local/share/orpara/matrix";
   }
}

// no validation, just set it, error will
// happen later
void MatrixScoreMethod::setDefaultPath(const string &path) {
   default_path = path;
}

/////////////// Constructors ///////////////
//

/*
 * Not making a copy of default protein matrix.
 * Waste of operation.
 *
MatrixScoreMethod::MatrixScoreMethod() 
   : path(default_path),
     name("blosum50"), matchS(5), mismatchS(-5),
     mat(), mtype(AMINO), 
     mins(numeric_limits<int>::max()), maxs(-1),
     words(0), wordSize(0)       
{ 
   // make a copy of the default matrix.
   alpha = 999999999;
   for (int i=0; i<26; i++) {
      for (int j=0; j<26; j++) {
         mat[i][j]=default_mat[i][j];
         if (mat[i][j] < mins) mins = mat[i][j];
         if (mat[i][j] > maxs) maxs = mat[i][j];
      }
   }
   setDefaultGapParameter();
}

MatrixScoreMethod::MatrixScoreMethod(const string& matrixName, bool nucMatrixScoreMethod) 
   : path(default_path), name(matrixName), 
           matchS(5), mismatchS(-5),
           mtype(NUCLEIC), 
           mins(numeric_limits<int>::max()), maxs(-1),
           words(0), wordSize(0) 
{ 
   read(); 
}

void MatrixScoreMethod::setDefaultGapParameter() {
   alpha = 2*getMinScore();
   beta = getMinScore()/3;
}

*/
//bool MatrixScoreMethod::read(bool nucMatrixScoreMethod) {
bool MatrixScoreMethod::read(int (*hashfunc)(char)) {
   string infile = path + "/" + str2upper(name);
   // I can use this part to decide which hash function to use
   // but for matrices not recognized by isProteinMatrix we are 
   // out of luck.
   if (isProteinMatrix(name)) {
#ifdef DEBUG
      //mtype = AMINO;
      cerr << "Reading amino acid matrix from file: " 
         << infile << " ...\n";
#endif
   }
   else if (isNucleicMatrix(name)) {
      //mtype = NUCLEIC;
#ifdef DEBUG
      cerr << "Reading nucleic acid matrix from file: "
        << infile << " ...\n";
#endif
   }
   else {
      if (name.empty()) {
         cerr << "Name is empty inside read!\n";
         throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Empty matrix name");
         //exit(1);
      }
      cerr << "matrix " << name << " is not known in the path: " << path << endl;
      //throw runtime_error("matrix: " + name + " unknown type in " + path);
      //exit(1);
      return false;
   }

   ifstream IN(infile.c_str());
   if (IN.fail()) {
      cerr << "Failed to open " << infile << " while trying to read scoring matrix\n";
      throw runtime_error("bad matrix file: " + infile);
      //exit(1);
   }
   string line;
   // read the comment and discard them
   getline(IN, line);
   while (!IN.eof() && line[0] == '#' || line.length()<2) {
      getline(IN, line);
   }
   unsigned int i;
   //char symb[27]; // all possible symbols for both AA and Base
   // read the first header line of the matrix
   vector<string> row = dissect(line, " \t");
   for (i=0; i<row.size(); i++) {
      if (row[i].size()>1) {
         cerr << "Violated the one-letter length" << row[i] << endl;
         //exit(1);
         return false;
      }
      symb[i]=row[i][0];
   }
   symb[i]='\0'; // residues symbols, total alphabet
   numsymbol=row.size();
   // for debug
   //cerr << numsymbol << " matrix letters\n"
   //   << symb << endl;

   getline(IN, line);
   // these two are initialized at construction.
   // For repeated reading they needs to be reset
   maxs=numeric_limits<int>::min();
   mins=numeric_limits<int>::max();
   int ss, r;
   // read the whole matrix
   while (!IN.eof() && !line.empty()) {
      row=dissect(line, " \t");
      if (row.size() != numsymbol + 1) {
         cerr << __FILE__ << ":" << __LINE__ << ":ERROR One row must have " << numsymbol << "elements\n";
         //exit(1);
         return false;
      }
      /*
      if (mtype == NUCLEIC) {
         r = hashbase(row[0][0]);
      }
      else if (mtype == AMINO) {
         r=aachar2num(row[0][0]);
      }
      else {
         cerr << "IDENTITY matrix has no real matrix data!\n";
         exit(1);
      }
      */
      r = hashfunc(row[0][0]);

      for (i=0; i<numsymbol; i++) {
         ss=atoi(row[i+1].c_str());
         if (ss > maxs) maxs=ss;
         if (ss < mins) mins=ss;
         /*
         if (mtype == NUCLEIC) 
            mat[r][hashbase(symb[i])]=ss; //atoi(tmp[i+1].c_str());
         else 
            mat[r][aachar2num(symb[i])]=ss; //atoi(tmp[i+1].c_str());
         */
         mat[r][hashfunc(symb[i])] = ss;
      }
      getline(IN, line);
   }
   //setDefaultGapParameter();
   return true;
}

void MatrixScoreMethod::deallocateWords() {
   if (words != nullptr) {
      for (size_t i=0; i<wordsArraySize; i++) {
         delete[] words[i];
      }
      delete[] words;
   }
}

MatrixScoreMethod::~MatrixScoreMethod() {
   deallocateWords();
}

void MatrixScoreMethod::setMinMaxScore() const {
   for (int i=0; i < getNumberOfAlphabet(); ++i) {
      for (int j=0; j < getNumberOfAlphabet(); ++j) {
         if (mat[i][j] < mins) mins = mat[i][j];
         if (mat[i][j] > maxs) maxs = mat[i][j];
      }
   }
}

int MatrixScoreMethod::getMinScore() const { 
   if (mins == numeric_limits<int>::max()) {
      setMinMaxScore();
   }
   return mins; 
}

int MatrixScoreMethod::getMaxScore() const { 
   if (maxs == -1) {
      setMinMaxScore();
   }
   return maxs; 
}

// copy constructor
MatrixScoreMethod::MatrixScoreMethod(const MatrixScoreMethod& mt) 
  : ScoreMethod(mt),
      path(mt.path), name(mt.name), 
      mat{0}, symb{' '}, numsymbol(mt.numsymbol),
      maxs(mt.maxs), mins(mt.mins), 
      words(nullptr), wordsArraySize(0), wordSize(0)
{
   setMatrix(mt.mat, 32);
   for (size_t i=0; i<32; ++i) {
      symb[i] = mt.symb[i];
   }
#ifdef DEBUG
   cerr << "called MatrixScoreMethod copy constructor\n";
#endif
}

MatrixScoreMethod::MatrixScoreMethod(MatrixScoreMethod&& mt)
   : ScoreMethod(std::move(mt)),
     path(std::move(mt.path)), name(std::move(mt.name)),
     mat{0}, symb{' '}, numsymbol(mt.numsymbol), 
     mins(mt.mins), maxs(mt.maxs),
     words(nullptr), wordsArraySize(0), wordSize(0)
{
   setMatrix(mt.mat, 32);
   for (size_t i=0; i<32; ++i) {
      symb[i] = mt.symb[i];
   }
   mt.deallocateWords();
   mt.numsymbol=0;

}

MatrixScoreMethod& MatrixScoreMethod::operator=(const MatrixScoreMethod& mt) {
   if (this != &mt) {
      ScoreMethod::operator=(mt);
      path=mt.path;
      name=mt.name;
      for (size_t i=0; i<32; ++i) {
         symb[i] = mt.symb[i];
      }
      mins=mt.mins;
      maxs=mt.maxs;
      numsymbol=mt.numsymbol;
      words=nullptr;
      wordsArraySize=0;
      wordSize=0;
      setMatrix(mt.mat, 32);
   }
   return *this;
}

MatrixScoreMethod& MatrixScoreMethod::operator=(MatrixScoreMethod&& mt) {
   if (this != &mt) {
      ScoreMethod::operator=(std::move(mt));
      path=std::move(mt.path);
      name=std::move(mt.name);
      for (size_t i=0; i<32; ++i) {
         symb[i] = mt.symb[i];
      }
      mins=mt.mins;
      maxs=mt.maxs;
      numsymbol=mt.numsymbol;
      words=nullptr;
      wordsArraySize=0;
      wordSize=0;
      setMatrix(mt.mat, 32);
      mt.deallocateWords();
   }
   return *this;
}

void MatrixScoreMethod::show(ostream &ous, int (*hashfunc)(char)) const {
   ScoreMethod::show(ous);
   ous << "Name: " << name << "\nPath: " << path << endl;
   int i;
   ous << "   ";
   int maxh=-1;
   char hashedsymb[32];
   for (i=0; i<32; ++i) hashedsymb[i] = ' ';
   for (i=0; i<numsymbol; i++) {
      ous << setw(2) << symb[i] << " ";
      hashedsymb[hashfunc(symb[i])] = symb[i];
      if (hashfunc(symb[i]) > maxh) {
         maxh = hashfunc(symb[i]);
      }
   } 
   ous << endl << "hashed position" << endl << "   ";
   for (i=0; i <= maxh; i++) {
      ous << setw(2) << hashedsymb[i] << " ";
   } 
   ous << endl;
   for (i=0; i <= maxh; i++) {
      ous << hashedsymb[i] << " {";
      for (int j=0; j <= maxh; j++) {
         ous << setw(2) << mat[i][j] << " ";
      }
      ous << "}\n";
   }
}

void MatrixScoreMethod::growWord(const vector<string> &in, vector<string> &ou) const {
   if (!ou.empty()) ou.clear();
   //cerr << "growing words from " << in.size() << endl;
   for (int i=0; i<in.size(); i++) {
      for (int j=0; j<numsymbol; j++) {
         ou.push_back(in[i] + symb[j]);
     }
   }
}

void MatrixScoreMethod::expandWord(char** &in, int &insize, int &ws) const {
   int ousize=insize*numsymbol;
   int i, k=0;
   char** ou=new char*[ousize];
   for (i=0; i<insize; i++) {
      for (int j=0; j<numsymbol; j++) {
         ou[k]=new char[wordSize+1];
         memcpy(ou[k], in[i], ws);
         ou[k][ws]=symb[j];
         ou[k][ws+1]='\0';
         ++k;
     }
     delete[] in[i];
   }
   delete[] in;
   ++ws;
   in=ou;
   insize=ousize;
}

// non-member function
void expandCombination(char** &buff, int &buffsize, int &ws, const char* alphabet, int alphabetsize) {
   int ousize=buffsize*alphabetsize;
   int i, k=0;
   char** ou=new char*[ousize];
   for (i=0; i<buffsize; i++) {
      for (int j=0; j<alphabetsize; j++) {
         ou[k]=new char[ws+2];
         memcpy(ou[k], buff[i], ws);
         ou[k][ws]=alphabet[j];
         ou[k][ws+1]='\0';
         ++k;
     }
     delete[] buff[i];
   }
   delete[] buff;
   ++ws;
   buff=ou;
   buffsize=ousize;
}

// failed test
char** MatrixScoreMethod::allwords(int ws) const {
   int i;
   // had what we wanted
   if (ws == wordSize) return words;
   // discard old memory
   if (words != nullptr) {
      for (i=0; i<wordsArraySize; i++) {
         delete[] words[i];
      }
      delete[] words;
   }
   cerr << "Making all words of word size: " << ws << endl;
   char AA[32];
   char *p=AA;
   int numAA=0;
   for (i=0; i<numsymbol; i++) {
      if (symb[i] != 'X' && symb[i] != '*' && symb[i] != 'Z'
            && symb[i] != 'B' && symb[i] != 'U') {
         *p=symb[i];
         ++p;
         ++numAA;
      }
   }
   *p='\0';

   words=new char*[numAA];
   for (i=0; i<numAA; i++) {
      words[i]=new char[2];
      words[i][0] = AA[i];
      words[i][1] = '\0';
   }
   wordsArraySize=numAA;
   int currws=1;
   for (i=0; i<ws-1; i++) {
      //expandWord(words, wordsArraySize, currws);
      expandCombination(words, wordsArraySize, currws, AA, numAA);
   }
   wordSize=ws;
   return words;
}

void MatrixScoreMethod::showWords() const {
   if (words == nullptr) {
      cerr << "words is empty!\n";
      return;
   }
   for (int i=0; i<wordsArraySize; i++) {
      cout << words[i] << " ";
      if ((i+1)%numsymbol == 0) cout << endl;
   }
   cout << endl;
}
      

void MatrixScoreMethod::getWords(vector<string> &words, int ws) const {
   vector<string> seed;
   int i;
   for (i=0; i<numsymbol; i++) {
      seed.push_back(string(1, symb[i]));
   }
   //cerr << ws << " word size inside getWords()\n";
   for (i=0; i<ws-1; i++) {
      //cerr << i << " th turn\n";
      growWord(seed, words);
      seed=words;
   }
}

int MatrixScoreMethod::similarWord_debug(ostream &ous, const char* w, int ws, int cutoff, float fractioncutoff) const {
   if (words == 0) {
      allwords(ws);
   }

   if (ws != wordSize) {
      cerr << "word must be the same size: " << wordSize << endl;
      exit(1);
   }
   //cout << "Words similar to " << w << " cutoff=" << cutoff << endl;
   int maxscore=0, selfscore;
   vector<pair<char*, int> > wordscore; 
   int i,j,k, score;
   for (k=0; k<wordsArraySize; k++) {
      score=0;
      for (i=0; i<wordSize; i++) score += lookup(w[i], words[k][i]);
      if (strcmp(w, words[k])) {
         if (score > maxscore) maxscore=score;
         if (score > cutoff)
            wordscore.push_back(make_pair(words[k], score));
      }
      else selfscore=score;
   }
   vector<pair<char*, int> > neighborwords; 
   for (int j=0; j<wordscore.size(); j++) {
      if (wordscore[j].second > fractioncutoff*selfscore) {
         neighborwords.push_back(wordscore[j]);
      }
   }
   if (!neighborwords.empty()) {
      ous << "Words similar to " << w << " cutoff=" << cutoff << endl;
      ous << "selfscore=" << selfscore << " maxscore=" << maxscore << endl;
      for (i=0; i<neighborwords.size(); i++) {
         ous << neighborwords[i].first << "->" 
            << neighborwords[i].second << "; ";
      }
      ous << "\n------------------------------------\n";
   }
   return neighborwords.size();
}

int MatrixScoreMethod::similarWord(vector<char*> &neighbor, const char* w, int ws, float cutoff, float fractioncutoff) const {
   if (words == 0) {
      allwords(ws);
   }
   if (ws != wordSize) {
      cerr << "word must be the same size: " << wordSize << endl;
      exit(1);
   }
   //int maxscore=0, selfscore;
   int selfscore;
   vector<pair<char*, int> > wordscore; 
   int i,j,k, score;
   for (k=0; k<wordsArraySize; k++) {
      score=0;
      for (i=0; i<wordSize; i++) score += lookup(w[i], words[k][i]);
      if (strcmp(w, words[k])) {
         //if (score > maxscore) maxscore=score;
         if (score > cutoff)
            wordscore.push_back(make_pair(words[k], score));
      }
      else selfscore=score;
   }
   if (!neighbor.empty()) neighbor.clear();
   //vector<pair<char*, int> > neighborwords; 
   for (int j=0; j<wordscore.size(); j++) {
      if (wordscore[j].second > fractioncutoff*selfscore) {
         neighbor.push_back(wordscore[j].first);
      }
   }
   /*
   if (!neighborwords.empty()) {
      ous << "Words similar to " << w << " cutoff=" << cutoff << endl;
      ous << "selfscore=" << selfscore << " maxscore=" << maxscore << endl;
      for (i=0; i<neighborwords.size(); i++) {
         ous << neighborwords[i].first << "->" 
            << neighborwords[i].second << "; ";
      }
      ous << "\n------------------------------------\n";
   }
   */
   return neighbor.size();
}

int MatrixScoreMethod::score(const int* x, const int* y, const int len) const {
   int s=0;
   for (int i=0; i<len; i++) {
      s += lookup(x[i], y[i]);
   }
   return s;
}

void MatrixScoreMethod::setMatrix(const int source[][32], const int size) {
   for (size_t i=0; i < size; i++) {
      for (size_t j=0; j < size; j++) {
         mat[i][j]=source[i][j];
      }
   }
}

void MatrixScoreMethod::copyFrom(const string &matName, 
      const char matSymb[32], const int matSize, const int scoreMat[][32]) {
   setName(matName);
   numsymbol = matSize;
   setAlphabet(matSymb);
   setMatrix(scoreMat, matSize);
}

/////////// ProteinScoringMethod //////////////////////////

// silly enough, the following is the same as = "ARND..."!
// I use this one for clarity, This can be used as a C-string because of the
// terminator. 24 alphabets, one terminator, it should be 25 for copying
// purpose.
string ProteinScoreMethod::default_name = string("blosum50");
char ProteinScoreMethod::default_symb[32] = {'A', 'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',  'B',  'Z',  'X',  '*', '\0'};

int ProteinScoreMethod::default_numsymbol = 24;
// This is a default protein matrix.
int ProteinScoreMethod::default_mat[32][32] = {
   /* this is the default blosum50 matrix in computer format*/
{5, -2, -1, -2, -1, -3, 0, -2, -1, 0, -1, -2, -1, -1, 0, -1, -1, -2, 1, 0, 0, 0, -3, -1, -2, -1, -5, 0, 0, 0 , 0, 0, },
{-2, 5, -3, 5, 1, -4, -1, 0, -4, 0, 0, -4, -3, 4, 0, -2, 0, -1, 0, 0, 0, -4, -5, -1, -3, 2, -5, 0, 0, 0, 0, 0, },
{-1, -3, 13, -4, -3, -2, -3, -3, -2, 0, -3, -2, -2, -2, 0, -4, -3, -4, -1, -1, 0, -1, -5, -2, -3, -3, -5, 0,  0, 0, 0, 0, },
{-2, 5, -4, 8, 2, -5, -1, -1, -4, 0, -1, -4, -4, 2, 0, -1, 0, -2, 0, -1, 0, -4, -5, -1, -3, 1, -5, 0, 0, 0, 0, 0, },
{-1, 1, -3, 2, 6, -3, -3, 0, -4, 0, 1, -3, -2, 0, 0, -1, 2, 0, -1, -1, 0, -3, -3, -1, -2, 5, -5, 0, 0, 0, 0,  0, },
{-3, -4, -2, -5, -3, 8, -4, -1, 0, 0, -4, 1, 0, -4, 0, -4, -4, -3, -3, -2, 0, -1, 1, -2, 4, -4, -5, 0, 0, 0,  0, 0, },
{0, -1, -3, -1, -3, -4, 8, -2, -4, 0, -2, -4, -3, 0, 0, -2, -2, -3, 0, -2, 0, -4, -3, -2, -3, -2, -5, 0, 0, 0, 0, 0, },
{-2, 0, -3, -1, 0, -1, -2, 10, -4, 0, 0, -3, -1, 1, 0, -2, 1, 0, -1, -2, 0, -4, -3, -1, 2, 0, -5, 0, 0, 0, 0 , 0, },
{-1, -4, -2, -4, -4, 0, -4, -4, 5, 0, -3, 2, 2, -3, 0, -3, -3, -4, -3, -1, 0, 4, -3, -1, -1, -3, -5, 0, 0, 0 , 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{-1, 0, -3, -1, 1, -4, -2, 0, -3, 0, 6, -3, -2, 0, 0, -1, 2, 3, 0, -1, 0, -3, -3, -1, -2, 1, -5, 0, 0, 0, 0,  0, },
{-2, -4, -2, -4, -3, 1, -4, -3, 2, 0, -3, 5, 3, -4, 0, -4, -2, -3, -3, -1, 0, 1, -2, -1, -1, -3, -5, 0, 0, 0 , 0, 0, },
{-1, -3, -2, -4, -2, 0, -3, -1, 2, 0, -2, 3, 7, -2, 0, -3, 0, -2, -2, -1, 0, 1, -1, -1, 0, -1, -5, 0, 0, 0, 0, 0, },
{-1, 4, -2, 2, 0, -4, 0, 1, -3, 0, 0, -4, -2, 7, 0, -2, 0, -1, 1, 0, 0, -3, -4, -1, -2, 0, -5, 0, 0, 0, 0, 0 , },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{-1, -2, -4, -1, -1, -4, -2, -2, -3, 0, -1, -4, -3, -2, 0, 10, -1, -3, -1, -1, 0, -3, -4, -2, -3, -1, -5, 0,  0, 0, 0, 0, },
{-1, 0, -3, 0, 2, -4, -2, 1, -3, 0, 2, -2, 0, 0, 0, -1, 7, 1, 0, -1, 0, -3, -1, -1, -1, 4, -5, 0, 0, 0, 0, 0 , },
{-2, -1, -4, -2, 0, -3, -3, 0, -4, 0, 3, -3, -2, -1, 0, -3, 1, 7, -1, -1, 0, -3, -3, -1, -1, 0, -5, 0, 0, 0,  0, 0, },
{1, 0, -1, 0, -1, -3, 0, -1, -3, 0, 0, -3, -2, 1, 0, -1, 0, -1, 5, 2, 0, -2, -4, -1, -2, 0, -5, 0, 0, 0, 0, 0, },
{0, 0, -1, -1, -1, -2, -2, -2, -1, 0, -1, -1, -1, 0, 0, -1, -1, -1, 2, 5, 0, 0, -3, 0, -2, -1, -5, 0, 0, 0, 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{0, -4, -1, -4, -3, -1, -4, -4, 4, 0, -3, 1, 1, -3, 0, -3, -3, -3, -2, 0, 0, 5, -3, -1, -1, -3, -5, 0, 0, 0,  0, 0, },
{-3, -5, -5, -5, -3, 1, -3, -3, -3, 0, -3, -2, -1, -4, 0, -4, -1, -3, -4, -3, 0, -3, 15, -3, 2, -2, -5, 0, 0 , 0, 0, 0, },
{-1, -1, -2, -1, -1, -2, -2, -1, -1, 0, -1, -1, -1, -1, 0, -2, -1, -1, -1, 0, 0, -1, -3, -1, -1, -1, -5, 0, 0, 0, 0, 0, },
{-2, -3, -3, -3, -2, 4, -3, 2, -1, 0, -2, -1, 0, -2, 0, -3, -1, -1, -2, -2, 0, -1, 2, -1, 8, -2, -5, 0, 0, 0 , 0, 0, },
{-1, 2, -3, 1, 5, -4, -2, 0, -3, 0, 1, -3, -1, 0, 0, -1, 4, 0, 0, -1, 0, -3, -2, -1, -2, 5, -5, 0, 0, 0, 0, 0, },
{-5, -5, -5, -5, -5, -5, -5, -5, -5, 0, -5, -5, -5, -5, 0, -5, -5, -5, -5, -5, 0, -5, -5, -5, -5, -5, 1, 0, 0, 0, 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }

};

ProteinScoreMethod::ProteinScoreMethod(const string& matrixName)
   : MatrixScoreMethod() 
{ 
   if (!isProteinMatrix(matrixName)) 
      cerr << "WARN: " << matrixName << " is not known protein matrix!\n"; 
   setName(matrixName); 
   MatrixScoreMethod::read(&aachar2num); 
}  

bool ProteinScoreMethod::read(const string &p, const string &n) {
   setPath(p);
   setName(n);
   return MatrixScoreMethod::read(&aachar2num);
}

ProteinScoreMethod::ProteinScoreMethod(const string& matrixName, int go, int ge)
   : MatrixScoreMethod(go, ge) 
{ 
   if (!isProteinMatrix(matrixName))
      cerr << "WARN: " << matrixName << " is not known protein matrix\n";
   setName(matrixName); 
   MatrixScoreMethod::read(&aachar2num); 
}

void ProteinScoreMethod::use(const string &matName) { 
   if (!isProteinMatrix(matName))
      cerr << "WARN: " << matName << " is not known protein matrix\n";
   setName(matName); 
   MatrixScoreMethod::read(&aachar2num); 
}


////////////// NucleicScoreMethod class ///////////////

string NucleicScoreMethod::default_name = string("NUC.4.4");
char NucleicScoreMethod::default_symb[32] = "ATGCSWRYKMBVHDN";
int NucleicScoreMethod::default_numsymbol = 15;
int NucleicScoreMethod::default_mat[32][32] = {
   //A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
   { 5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2},
   {-4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2},
   {-4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2},
   {-4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2},
   {-4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1},
   { 1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1},
   { 1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1},
   {-4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1},
   {-4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1},
   { 1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1},
   {-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1},
   {-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1},
   {-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1},
   {-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1},
   {-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

bool NucleicScoreMethod::read(const string &p, const string &n) {
   setPath(p);
   setName(n);
   return MatrixScoreMethod::read(&hashbase);
}

NucleicScoreMethod::NucleicScoreMethod(const string& matrixName)
   : MatrixScoreMethod(matrixName) 
{ 
   if (!isNucleicMatrix(matrixName)) {
      //cerr << "WARN: " << matrixName << " is not known Nucleic matrix!\n";
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR: " + matrixName + " unknown nucleic");
   }
   MatrixScoreMethod::read(&hashbase); 
}

NucleicScoreMethod::NucleicScoreMethod(const string& dir, const string& matrixName)
   : MatrixScoreMethod(dir, matrixName) 
{ 
   if (!isNucleicMatrix(matrixName)) {
      //cerr << "WARN: " << matrixName << " is not known Nucleic matrix!\n";
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR: " + matrixName + " unknown nucleic");
   }
   MatrixScoreMethod::read(&hashbase); 
}

NucleicScoreMethod::NucleicScoreMethod(const string& matrixName, int go, int ge)
   : MatrixScoreMethod(go, ge) 
{ 
   if (!isNucleicMatrix(matrixName)) {
      cerr << "WARN: " << matrixName << " is not known Nucleic matrix!\n";
      exit(1);
   }
   setName(matrixName); 
   MatrixScoreMethod::read(&hashbase); 
}

void NucleicScoreMethod::use(const string &matName) { 
   if (!isNucleicMatrix(matName)) {
      cerr << "WARN: " << matName << " is not known Nucleic matrix!\n";
      exit(1);
   }
   setName(matName); 
   MatrixScoreMethod::read(&hashbase); 
}
}
