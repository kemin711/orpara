#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include "bioseq.h"

using namespace std;

namespace orpara {
//int aachar2num(char a);
//declared in bioseq.h
/*
int max(int i1, int i2, int i3) {
   return max(max(i1,i2),i3);
}
int max(int i1, int i2) {
   if (i1>i2) return i1;
   return i2;
}
*/


void expandCombination(char** &buff, int &buffsize, int &ws, const char* alphabet, int alphabetsize);

enum matrixType { NUCLEIC = 0, AMINO = 1, IDENTITY = 2 };

/**
 * Once the Matrix is constructed, it should be passed around as pointers.
 *
 * Not using singleton pattern because you may need several different
 * matrices in one application.  The programmer has to make sure not to
 * instanciate too many of the same kind.
 *
 * There are still a lots of problems to fix in this class.
 * This class is replaced by scorematrix.h ScoreMethod class hiarche.
 */
class Matrix {
   public:
      /** 
       * Default matrix constructor.
       *
       * The default matrix is for protein.
       * Makes a copy of the default matrix in memory.
       */
      Matrix();
      /*: path(default_path),
            name("blosum50"), matchS(5), mismatchS(-5),
            mat(), mtype(AMINO), words(0), wordSize(0)       
      { 
         // make a copy of the default matrix.
         for (int i=0; i<26; i++) for (int j=0; j<26; j++) mat[i][j]=default_mat[i][j];
      }
      */
      /** 
       * @param matrixName a known scoring matrix name
       * @param nucMatrix whether it is a matrix for nucleotide default is not.
       */
      Matrix(const string& matrixName, bool nucMatrix=false);
      /*
         : path(default_path), name(matrixName), 
           matchS(5), mismatchS(-5),
           mtype(NUCLEIC), words(0), wordSize(0) 
      { read(); }
      */
      /**
       * copy constructor
       */
      Matrix(const Matrix& mt);

      /** construct a matrix that has only two parameters */
      Matrix(int m, int mis) : matchS(m), mismatchS(mis), mtype(IDENTITY), alpha(2*mis), beta(mis/3),
         mins(mis), maxs(m) { }

      ~Matrix();

      /**
       * Get the minimum score of the whole matrix.
       */
      int getMinScore() const;
      /**
       * @return Maximum score of the whole matrix.
       */
      int getMaxScore() const;
      /**
       * Normally we don't copy matrices, because we usually
       * only need one for the whole application!
       * We usaully pass pointers around.
       */
      Matrix& operator=(const Matrix& mt);
      /** Debug function
       *  to show the matrix in the actual array format
       * output to cout
       */
      void show() const;

      /**
       * returns the alphabet of the all amino acid symbols
       * in the matrix
       */
      const char* getAlphabet() const { return aas; }
      /** returns the number of symbols in the matrix
       * for protein matrix this is usually 22, but with
       * selenocystein, this could be 23. 
       *
       * For nucleotide acids, this should be 5 with N
       * being the 5th symbol.
       */
      int getNumberOfAlphabet() const { return numsymbol; }

      /** 
       * put all possible words of length ws into the vector words
       */
      void getWords(vector<string> &words, int ws) const;
      void growWord(const vector<string> &in, vector<string> &ou) const;
      void expandWord(char** &in, int &insize, int &ws) const;
      /**
       * return a pointer to all the words generated in all combination.
       * Symbol B Z X * U should not be used to build words.
       */
      char** allwords(int ws) const;
      int getNumberOfWords() const { return wordsArraySize; }
      int getCurrentWordSize() const { return wordSize; }
      void showWords() const; // debug function
      /**
       * debug show all similar words and scores
       * input w  The input word, 
       * this function finds all similar words and write them to 
       * the output stream.
       *
       * For the real version, some words never have neighboring words,
       * We can use a set, pre-load the words into a set, for 
       * speedy look up.  This is for the future.
       */
      int similarWord_debug(ostream &ous, const char* w, int ws, int cutoff=10, float fractioncutoff=0.75) const;
      /**
       * neighbor will be cleared first if it is not empty.
       */
      int similarWord(vector<char*> &neighbor, const char* w, int ws, float cutoff=11, float fractioncutoff=0.8) const;

      /**
       * @return true for success, false for failure.
       * 
       * Type of matrix is deduced automatically.
       * NUCLEIC use hashbase() function defined in codon.h to
       * convert Base Char to integer index position.
       * AMINO use aachar2num to index position, other than the generic 
       *     char X-'A' function
       * */
      bool read();
      // no longer need to tell what type of matrix. Done automatically.
      //bool read(bool nucMatrix=false);
      /**
       * For using any user defined matrix.
       * Matric name must exists in the class level variable.
       *  proteinMatrices[] or nucleicMatrices[];
       */
      bool read(const string &p, const string &n) { setPath(p); name=n; return read(); }

      /** read a new matrix from file. Expensive operation.
       */
      void setMatrix(const string &n) { name=n; read(); }
      /**
       * set nucleic acid matrix.
       * Deprecated.
       */
      void setMatrixNuc(const string &n) { name=n; mtype=NUCLEIC; read(); }
      /**
       * alpha = 2*minScore
       * beta = 0.5*alpah
       */
      void setDefaultGapParameter();
      void setMinMaxScore() const;

      void setPath(const string &p) { path=p; }
      int getMatchScore(const char c1, const char c2) {
         if (c1==c2) return matchS; return mismatchS; }
      /**
       * More efficient using code that is the result of the
       * getcode() method of bioseq.
       *
       * In debugged code, the bound checking shoould be removed.
       * */
      int lookup(int row, int col) const { 
         if (mtype == IDENTITY) {
            if (row == col) return matchS;
            else return mismatchS;
         }
         /*
         if (row > 26 || col > 26) {
            cerr << "Matrix index out of bound " << row << ", " << col << endl;
            exit(1);
         }
         */
         return mat[row][col]; 
      }

      /**
       * use the aachar2num function in bioseq to
       * convert letters to numbers with A as zero, B as 1, ...
       * All letters must be upper case
       * Depends on the type of matrix, if nucleic acid matrix, then
       * use hashbase function,
       * else use aachar2num function.
       */
      int lookup(char row, char col) const { 
         if (mtype == NUCLEIC) {
            return mat[hashbase(row)][hashbase(col)]; 
         }
         else if (mtype == AMINO) {
            return mat[aachar2num(row)][aachar2num(col)]; 
         }
         else {
            if (row == col) return matchS;
            else return mismatchS;
         }
      }

      /**
       * give two array of int code for sequence, it return
       * the score for this segment of sequence alignment
       */
      int score(const int* x, const int* y, const int len) const;
      bool isNucleicAcid() const { return mtype == NUCLEIC; }
      bool isAminoAcid() const { return mtype == AMINO; }
      /**
       * Will not enforce both parameters being negative.
       */
      void setGapParameters(const int go, const int ge) {
         alpha=go; beta=ge; }
      void setGapInsert(const int go) { alpha = go; }
      void setGapExtend(const int ge) { beta = ge; }
      /**
       * @return the gap opening score (alpha).
       * used in gap(g) = alpha + (g-1)*beta function.
       */
      int getGapInsert() const { return alpha; }
      int getGapExtend() const { return beta; }

      /**
       * set the defalt maxtrix path to the matrix files
       * This method will set the default matrix path to 
       * $HOME/src/proj/seqaln/matrix
       * Without checking for errors.
       */
      static void setDefaultPath();
      static void setDefaultPath(const string &path);
      static bool isProteinMatrix(const string &matrixName);
      static bool isNucleicMatrix(const string &matrixName);

   private:
      /**
       * The path to the matrix directory
       */
      string path;
      /**
       * Name of the matrix
       */
      string name; // Blosum50, Blosum62, etc
      /** match score for identity matrix.
       * Identity matrix is very useful for dealing with sequence assemblies
       * where the symbol has no biological meaning.
       * Alternatively, these two values take the first two spaces of the 
       * mat store, thus eliminating the extra variable.
       * */
      int matchS;
      /** mismatch score for identity matrix */
      int mismatchS;
      /**
       * The core data structure for storing the actual information.
       */
      int mat[32][32]; // use only 27 of the elements
      //bool isNuc;
      /**
       * matrix type should be determined by this class
       * according to some static map.
       */
      matrixType mtype;
      /**
       * the gap score, this should have something to do 
       * with the matrix
       * g is the length of the gap.
       * gap(g)=alpha + (g-1)*beta
       * The gap function.
       * this is used in the recursion.
       * Both alpha and beta should be negative values.
       */
      int alpha, beta;
      /** 
       * minimum and maximum score in the matrix
       * This is a cache value to prevent repeated
       * calculation from the matrix.
       * maxs start with -1.
       * mins start with max_int
       */
      mutable int mins, maxs;

      /**
       * This is a default protein matrix initialized in the implementation.
       */
      static int default_mat[32][32];
      static string default_path;
      /**
       * names of protein matrices in lower case
       */
      static const char* proteinMatrices[];
      /**
       * names of nucleic acid matrices in lower case
       */
      static const char* nucleicMatrices[];

      /** residue symbols in temporary working space 
       * This will apply to both nucleic acid and amino acid.
       * For identical matrix this should be undefined.
       * The actual elements stored is defined by 
       * numsymbol value.
       * */
      char aas[32]; // amino acid symbols, at most 26 so we had enough
      /** 
       * size of aas array, go with aas.
       */
      int numsymbol; // the size of aas array

      mutable char **words;
      mutable int wordsArraySize;
      mutable int wordSize;  // recoreds the size of words
};
}
#endif
