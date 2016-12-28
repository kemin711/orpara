#include "matrix.h"
#include <vector>
#include "strformat.h"
#include <cstring>
#include <iomanip>

//#include <cstdlib>

namespace orpara {
// helper functioion
// moved to bioseq.cpp
//int aachar2num(char a) {
 //  int i=a-'A';
  // if (i < 0 ) return 26;
//   return i;
//}
string Matrix::default_path="/home/zhouke/src/proj/seqaln/matrix";
const char* Matrix::proteinMatrices[]  = {
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

const char* Matrix::nucleicMatrices[] = {
   "DNA", "DNAelem", "NUC.4.4" 
};

bool Matrix::isProteinMatrix(const string &matrixName) {
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

bool Matrix::isNucleicMatrix(const string &matrixName) {
   for (int i = 0; i<3; ++i) {
      if (nucleicMatrices[i] == matrixName) return true;
   }
   return false;
}

void Matrix::setDefaultPath() {
   string pathToMatrix(getenv("HOME"));
   pathToMatrix += "/src/proj/seqaln/matrix";
   default_path=pathToMatrix;
}
void Matrix::setDefaultPath(const string &path) {
   default_path = path;
}


   
// This is a default protein matrix.
int Matrix::default_mat[32][32] = {
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

Matrix::Matrix() 
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

Matrix::Matrix(const string& matrixName, bool nucMatrix) 
   : path(default_path), name(matrixName), 
           matchS(5), mismatchS(-5),
           mtype(NUCLEIC), 
           mins(numeric_limits<int>::max()), maxs(-1),
           words(0), wordSize(0) 
{ 
   read(); 
}

void Matrix::setDefaultGapParameter() {
   alpha = 2*getMinScore();
   beta = getMinScore()/3;
}

Matrix::~Matrix() {
   if (mtype == IDENTITY) return;
   if (words != 0) {
      for (int i=0; i<wordsArraySize; i++) {
         delete[] words[i];
      }
      delete[] words;
   }
}

void Matrix::setMinMaxScore() const {
   for (int i=0; i < getNumberOfAlphabet(); ++i) {
      for (int j=0; j < getNumberOfAlphabet(); ++j) {
         if (mat[i][j] < mins) mins = mat[i][j];
         if (mat[i][j] > maxs) maxs = mat[i][j];
      }
   }
}

int Matrix::getMinScore() const { 
   if (mins == numeric_limits<int>::max()) {
      setMinMaxScore();
   }
   return mins; 
}

int Matrix::getMaxScore() const { 
   if (maxs == -1) {
      setMinMaxScore();
   }
   return maxs; 
}

//bool Matrix::read(bool nucMatrix) {
bool Matrix::read() {
   string infile = path + "/" + name;
   if (isProteinMatrix(name)) {
      mtype = AMINO;
      //cerr << "Reading amino acid matrix from file: " 
      //   << infile << " ...\n";
   }
   else if (isNucleicMatrix(name)) {
      mtype = NUCLEIC;
      //cerr << "Reading nucleic acid matrix from file: "
      //   << infile << " ...\n";
   }
   else {
      cerr << "matrix " << name << " is now known in the path: " << path << endl;
      exit(1);
   }

   /* automatic assignment of matrix type, getting rid of 
    * unnecessary parameters
   if (nucMatrix) {
      mtype = NUCLEIC;
   }
   else mtype = AMINO;
   */

   // initialize to zero
   //int i,j;
   /* this is not needed because the default is 0
   for (i=0; i<28; i++) {
      for (j=0; j<28; j++) mat[i][j]=0;
   }
   */
   ifstream IN(infile.c_str());
   if (IN.fail()) {
      cerr << "Failed to open " << infile << endl;
      exit(1);
   }
   string line;
   // read the comment and discard them
   getline(IN, line);
   while (line[0] == '#' || line.length()<2) {
      getline(IN, line);
   }
   unsigned int i;
   //char aas[27]; // all possible symbols for both AA and Base
   // read the first header line of the matrix
   vector<string> row = dissect(line, " \t");
   for (i=0; i<row.size(); i++) {
      if (row[i].size()>1) {
         cerr << "Violated the one-letter length" << row[i] << endl;
         exit(1);
      }
      aas[i]=row[i][0];
   }
   aas[i]='\0'; // residues symbols, total alphabet
   numsymbol=row.size();
   // for debug
   //cerr << numsymbol << " matrix letters\n"
   //   << aas << endl;

   getline(IN, line);
   maxs=-99999999;
   mins=numeric_limits<int>::max();
   int ss, r;
   // read the whole matrix
   while (!IN.eof() && !line.empty()) {
      row=dissect(line, " \t");
      if (row.size() != numsymbol + 1) {
         cerr << "One row must have " << numsymbol << "elements\n";
         exit(1);
      }
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

      for (i=0; i<numsymbol; i++) {
         ss=atoi(row[i+1].c_str());
         if (ss > maxs) maxs=ss;
         if (ss < mins) mins=ss;
         if (mtype == NUCLEIC) 
            mat[r][hashbase(aas[i])]=ss; //atoi(tmp[i+1].c_str());
         else 
            mat[r][aachar2num(aas[i])]=ss; //atoi(tmp[i+1].c_str());
      }
      getline(IN, line);
   }
   setDefaultGapParameter();
}

Matrix::Matrix(const Matrix& mt) 
  : path(mt.path), name(mt.name), matchS(mt.matchS), mismatchS(mt.mismatchS),
      mtype(mt.mtype), alpha(mt.alpha), beta(mt.beta),
      maxs(mt.maxs), mins(mt.mins), numsymbol(mt.numsymbol),
      words(0), wordsArraySize(0), wordSize(0)
{
   for (int i=0; i<32; i++) {
      for (int j=0; j<32; j++) {
         mat[i][j]=mt.mat[i][j];
      }
   }
   strcpy(aas, mt.aas);
}

Matrix& Matrix::operator=(const Matrix& mt) {
   if (this != &mt) {
      path=mt.path;
      name=mt.name;
      matchS = mt.matchS;
      mismatchS = mt.mismatchS;
      for (int i=0; i<32; i++) {
         for (int j=0; j<32; j++) {
            mat[i][j]=mt.mat[i][j];
         }
      }
      alpha = mt.alpha;
      beta = mt.beta;
      maxs=mt.maxs;
      mins=mt.mins;
      strcpy(aas, mt.aas);
      numsymbol=mt.numsymbol;
      words=0;
      wordsArraySize=0;
      wordSize=0;
      mtype = mt.mtype;
   }
   return *this;
}

void Matrix::show() const {
   if (mtype == NUCLEIC) 
      cout << "Is nucleic acid matrix\n";
   else if (mtype == AMINO) 
      cout << "Is amino acid matrix\n";
   else {
      cout << "identity matrix.  Match: " << matchS
         << " mismatch: " << mismatchS << endl;
      return;
   }
   cout << "Name: " << name << "\nPath: " << path << endl;
   //for (int i=0; i<32; i++) {
   for (int i=0; i < numsymbol; i++) {
      cout << "{";
      for (int j=0; j < numsymbol; j++) {
         cout << setw(3) << mat[i][j] << ", ";
      }
      cout << "}\n";
   }
}

void Matrix::growWord(const vector<string> &in, vector<string> &ou) const {
   if (!ou.empty()) ou.clear();
   //cerr << "growing words from " << in.size() << endl;
   for (int i=0; i<in.size(); i++) {
      for (int j=0; j<numsymbol; j++) {
         ou.push_back(in[i] + aas[j]);
     }
   }
}

void Matrix::expandWord(char** &in, int &insize, int &ws) const {
   int ousize=insize*numsymbol;
   int i, k=0;
   char** ou=new char*[ousize];
   for (i=0; i<insize; i++) {
      for (int j=0; j<numsymbol; j++) {
         ou[k]=new char[wordSize+1];
         memcpy(ou[k], in[i], ws);
         ou[k][ws]=aas[j];
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

char** Matrix::allwords(int ws) const {
   int i;
   // had what we wanted
   if (ws == wordSize) return words;
   // discard old memory
   if (wordSize != 0) {
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
      if (aas[i] != 'X' && aas[i] != '*' && aas[i] != 'Z'
            && aas[i] != 'B' && aas[i] != 'U') {
         *p=aas[i];
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

void Matrix::showWords() const {
   for (int i=0; i<wordsArraySize; i++) {
      cout << words[i] << " ";
      if ((i+1)%numsymbol == 0) cout << endl;
   }
   cout << endl;
}
      

void Matrix::getWords(vector<string> &words, int ws) const {
   vector<string> seed;
   int i;
   for (i=0; i<numsymbol; i++) {
      seed.push_back(string(1, aas[i]));
   }
   //cerr << ws << " word size inside getWords()\n";
   for (i=0; i<ws-1; i++) {
      //cerr << i << " th turn\n";
      growWord(seed, words);
      seed=words;
   }
}

int Matrix::similarWord_debug(ostream &ous, const char* w, int ws, int cutoff, float fractioncutoff) const {
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

int Matrix::similarWord(vector<char*> &neighbor, const char* w, int ws, float cutoff, float fractioncutoff) const {
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

int Matrix::score(const int* x, const int* y, const int len) const {
   int s=0;
   for (int i=0; i<len; i++) {
      s += lookup(x[i], y[i]);
   }
   return s;
}
}
