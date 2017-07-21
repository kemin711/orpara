#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

// (c) 2012 Kemin Zhou at orpara.com

#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include "bioseq.h"

//#define DEBUG

using namespace std;

namespace orpara {
/**
 * Simplest score method and serve as the base class for all other
 * score methods include matrices for both protein and nucleic acids.
 * This is also the most effient class without any member variables.
 * Actually this class is very useful in dealing with sequences
 * for the purpose of identifying sequences without any regard for their
 * biological meaning. Next generation sequences should use this 
 * scoring method.
 *
 * Match score is 10 and mismatch score is -10.  These are hardcoded
 * and you cannot change. You can only change the gap open and extension
 * score.
 */
class ScoreMethod {
   public:
      /**
       * Default constructor.
       * default gap open -20, gap extend -1
       */
      ScoreMethod() : gapo(-20), gape(-1) { }
      ScoreMethod(int gopen, int gextend) : gapo(gopen), gape(gextend) { }
      /**
       * Copy constructor
       */
      ScoreMethod(const ScoreMethod &scm) : gapo(scm.gapo), gape(scm.gape) { }
      ScoreMethod(ScoreMethod &&scm) : gapo(scm.gapo), gape(scm.gape) { }
      ScoreMethod& operator=(const ScoreMethod &scm);
      ScoreMethod& operator=(ScoreMethod &&scm);
      /**
       * look up the score for two residues
       * @return the match score (10) if they are the same -10 if they are
       * different.
       */
      int lookup(const char r1, const char r2) const { if (r1 == r2) return 10; else return -10; }
      /**
       * look up the score for two interger representation of the two residues.
       * Same algorith as the char version.
       */
      int lookup(const int c1, const int c2) const { if (c1 == c2) return 10; else return -10; }

      int getMatch() const { return 10; }
      int getMismatch() const { return -10; }

      /**
       * @return -20 as the gap opening penalty score.
       */
      int getGapOpen() const { return gapo; }
      /**
       * Alias for getGapOpen()
       */
      int getGapInsert() const { return gapo; }
      /**
       * @return -1 as the gap extension penalty score.
       */
      int getGapExtend() const { return gape; }

      void setGapParameters(int gapopen, int gapextend) {
         gapo = gapopen; gape = gapextend; }
      void setGapOpen(int gopen) { gapo = gopen <= 0? gopen : -gopen; }
      void setGapExtend(int gextend) { gape = gextend <= 0? gextend : -gextend; }
      void show(ostream &ous) const;
      void show() const { show(cout); }

   protected:
      int gapo, gape;
};


/**
 * Identity score method. Same as the parent class except that the user of the
 * class can change the values of the match/mismatch parameters.  It simply
 * return one score for match, and another for mismatch. Gap penalty scores can
 * also be altered by the user, which is inherited from the parent.
 *
 */
class SimpleScoreMethod : public ScoreMethod {
   public:
      /**
       * Default constructor.
       * default match 10, mismatch -10.
       * The gap open is the default from the parent class -20,
       * gap extend -1.
       */
      SimpleScoreMethod() : ScoreMethod(), matchs(10), mismatchs(-10) { }
      SimpleScoreMethod(int match, int mismatch) 
         : ScoreMethod(2*mismatch, mismatch/8 <0? mismatch/8 : -1),
            matchs(match), mismatchs(mismatch) { }
      /**
       * complete constructor.
       * @param match match score
       * @param mismatch mismatch score
       * @param gapOpen score for opening a gap.
       * @param gapExtend score for extend a gap by one extra.
       */
      SimpleScoreMethod(int match, int mismatch, int gapOpen, int gapExtend) 
         : ScoreMethod(gapOpen, gapExtend), matchs(match), mismatchs(mismatch) { }
      /** copy constructor.
       */
      SimpleScoreMethod(const SimpleScoreMethod &meth) 
         : ScoreMethod(meth), matchs(meth.matchs), mismatchs(meth.mismatchs) { }
      /**
       * Assignment operator.
       */
      SimpleScoreMethod& operator=(const SimpleScoreMethod &mt);
      /**
       * look up the score for two residues
       * @return the match score (10) if they are the same -10 if they are
       * different.
       */
      int lookup(const char r1, const char r2) const { 
         if (r1 == r2) return matchs; else return mismatchs; }
      /**
       * look up the score for two interger representation of the two residues.
       * Same algorith as the char version.
       */
      int lookup(const int c1, const int c2) const { 
         if (c1 == c2) return matchs; else return mismatchs; }

      int getMatch() const { return matchs; }
      int getMismatch() const { return mismatchs; }
      void setMatch(int match) { matchs = match; }
      void setMismatch(int mis) { mismatchs = mis; }

   private:
      int matchs, mismatchs;
};

void expandCombination(char** &buff, int &buffsize, int &ws, const char* alphabet, int alphabetsize);

// we are not going to need this for performance
//
//enum matrixType { NUCLEIC = 0, AMINO = 1, IDENTITY = 2 };

/**
 * Once the Matrix is constructed, it should be passed around as pointers.
 *
 * Not using singleton pattern because you may need several different
 * matrices in one application.  The programmer has to make sure not to
 * instanciate too many of the same kind.
 *
 * This is more complicated based on a scoring table.
 * This class is protected and cannot be instantiated.
 * You must use the derive class: protein or nucleic acid.
 *
 * This class provide most of the functinalities for the derived classes.
 *
 * Typical usage is to set the default path of the matrix
 * then you can ignore the path part of the matrix.
 */
class MatrixScoreMethod : public ScoreMethod {
   public:
      /**
       * Set the defalt path to path to the scoring matrix files.
       * It will look for environmental variables SCORING_MATRIX_PATH first.
       * If found will use its value; otherwise will set default_path
       * to /usr/local/share/orpara/matrix without 
       * validation.
       */
      static void setDefaultPath();
      static const string& getDefaultPath() { return default_path; }
      /**
       * Set the default_path to path.
       * Use setPath for instance variable.
       */
      static void setDefaultPath(const string &path);
      static bool isProteinMatrix(const string &matrixName);
      /**
       * Helpter method to check whether this matrix
       * has been registered in this class.
       */
      static bool isNucleicMatrix(const string &matrixName);

   protected:
      /** 
       * Default matrix constructor.
       * Only intialize the gap paremeter according to the parent class
       * initializer. will set the default path to the matrix directory.
       */
      MatrixScoreMethod() 
         : ScoreMethod(), path(default_path), name(), mat{0}, 
           symb{' '}, numsymbol(0),
            mins(numeric_limits<int>::max()), maxs(-1),
            words(nullptr), wordsArraySize(0), wordSize(0) { }
      /**
       * Fixed gap open and gap extension.
       */
      MatrixScoreMethod(int go, int ge) 
         : ScoreMethod(go, ge), path(default_path), name(), mat{0}, 
           symb{' '}, numsymbol(0),
            mins(numeric_limits<int>::max()), maxs(-1),
            words(nullptr), wordsArraySize(0), wordSize(0) { }

      /** 
       * The caller should set the default matrix path
       * before calling this method. The matrix should
       * be installed in a typical directory. Then
       * the name correspnds to the files in the matrix
       * directory.
       * @param matrixName a known scoring matrix name which
       *   is also a matrix file name.
       * @param nucMatrixScoreMethod whether it is a matrix for nucleotide default is not.
       * Cannot initialize the matrix or read the actual matrix.
       * Require the type of matrix.
       */
      MatrixScoreMethod(const string& matrixName)
         : ScoreMethod(), path(default_path), name(matrixName), 
           mat{0}, symb{' '}, numsymbol(0), 
            mins(numeric_limits<int>::max()), maxs(-1),
            words(nullptr), wordsArraySize(0), wordSize(0) 
      { }

      /**
       * In case class level matrix path is not set,
       * you can use this method to create your new
       * matrix.
       *
       * @param dir path to matrix directory
       * @param matrixName name of the matrix.
       *
       * Note: This constructor Will not initialize the matrix use
       *   derived class to do it.
       */
      MatrixScoreMethod(const string& dir, const string& matrixName)
         : ScoreMethod(), path(dir), name(matrixName), 
           mat{0}, symb{' '}, numsymbol(0),
            mins(numeric_limits<int>::max()), maxs(-1),
            words(nullptr), wordsArraySize(0), wordSize(0) 
      { }

      MatrixScoreMethod(const string& matrixName, int go, int ge)
         : ScoreMethod(go, ge), path(default_path), name(matrixName), 
           mat{0}, symb{' '}, numsymbol(0),
            mins(numeric_limits<int>::max()), maxs(-1),
            words(nullptr), wordsArraySize(0), wordSize(0) 
      { }

   public:
      /**
       * copy constructor
       */
      MatrixScoreMethod(const MatrixScoreMethod& mt);
      /**
       * Move constructor
       */
      MatrixScoreMethod(MatrixScoreMethod&& mt);

      ~MatrixScoreMethod();

      void setName(const string &matName) { name = matName; }
      const string& getName() const { return name; }

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
      MatrixScoreMethod& operator=(const MatrixScoreMethod& mt);
      MatrixScoreMethod& operator=(MatrixScoreMethod&& mt);
      /** Debug function
       *  to show the matrix in the actual array format
       * output to cout
       */
      void show(int (*hashfunc)(char)) const { show(cout, hashfunc); }
      void show(ostream &ous, int (*hashfunc)(char)) const;

      /**
       * returns the alphabet of the all amino acid symbols
       * in the matrix
       */
      const char* getAlphabet() const { return symb; }
      /**
       * Set the alphabet used in this matrix.
       * The alphabet was also called symbols; I use the two terms
       * interchangeably.
       */
      void setAlphabet(const char symbols[32]) { strcpy(symb, symbols); }

      /** returns the number of symbols in the matrix
       * for protein matrix this is usually 22, but with
       * selenocystein, this could be 23. 
       *
       * For nucleotide acids, this should be 5 with N
       * being the 5th symbol.
       */
      int getNumberOfAlphabet() const { return numsymbol; }
      /**
       * set the internal numsymbol to na.
       */
      void setNumberOfAlphabet(const int na) { numsymbol = na; }

      /** Utility function.
       */
      void copyFrom(const string &matName, const char matSymb[32], 
            const int matSize, const int scoreMat[][32]);

      /** 
       * put all possible words of length ws into the vector words
       */
      void getWords(vector<string> &words, int ws) const;
      void growWord(const vector<string> &in, vector<string> &ou) const;
      void expandWord(char** &in, int &insize, int &ws) const;
      /**
       * return a pointer to all the words generated in all combination.
       * Symbol B Z X * U should not be used to build words.
       * TODO: this one failed test need to debug.
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

      // no longer need to tell what type of matrix. Done automatically.
      //bool read(bool nucMatrixScoreMethod=false);


      /**
       * set the path where the matrix text file is stored.
       */
      void setPath(const string &p) { path=p; }

      /**
       * More efficient using code that is the result of the
       * getcode() method of bioseq.
       *
       * In debugged code, the bound checking shoould be removed.
       * */
      int lookup(int row, int col) const { 
#ifdef DEBUG
         // this is for debug
         if (row > 26 || col > 26) {
            cerr << "MatrixScoreMethod index out of bound " << row << ", " << col << endl;
            exit(1);
         }
#endif
         return mat[row][col]; 
      }

      /**
       * give two array of int code for sequence, it return
       * the score for this segment of sequence alignment
       */
      int score(const int* x, const int* y, const int len) const;

   protected:

      /**
       * @return true for success, false for failure.
       * 
       * Type of matrix is deduced automatically.
       * NUCLEIC use hashbase() function defined in codon.h to
       * convert Base Char to integer index position.
       * AMINO use aachar2num to index position, other than the generic 
       *     char X-'A' function
       *
       * This method dpends on the type of symbol, and use
       * different hash functions. So it should be parameterized
       * to be used by different type of matrices.
       *
       * @param hashfunc is the pointer to the function hashbase or aachar2num
       * */
      bool read(int (*hashfunc)(char));
      /**
       * For using any user defined matrix.
       * Matric name must exists in the class level variable.
       *  proteinMatrices[] or nucleicMatrices[];
       */
      //bool read(const string &p, const string &n) { setPath(p); name=n; return read(); }

      /**
       * Copy from a source matrix into the mat filed of this object.
       * Note this method does not know the exact number of symbols
       * and assumes that the matrix to be copied has the same
       * number of alphabets as that of this one.
       * @param size is the size of the alphabet. For protein it is 24.
       *    Should not exceed 32.
       * @param source the source matrix that is allocated [32]x[32] 2-d int
       * array.
       */
      void setMatrix(const int source[][32], const int size=32);
      /**
       * Deallocate the words 2-dimensional array.
       */
      void deallocateWords();

      /**
       * The path to the matrix directory, this could be
       * an environment variable or given by user.
       */
      string path;
      /**
       * Name of the matrix
       *  such as Blosum50, Blosum62, etc
       */
      string name; // Blosum50, Blosum62, etc
      /**
       * The core data structure for storing the actual information.
       * Use 27 of the elements. 32 is a power of 2.
       */
      int mat[32][32]; 
      /** 
       * Alphabet
       * residue symbols in temporary working space 
       * This will apply to both nucleic acid and amino acid.
       * The actual elements stored is defined by 
       * numsymbol value.
       * */
      char symb[32]; 
      /** 
       * size of symb array, go with symb.
       * But the size of matrix could be bigger because amino acids
       * has no symbol for B, J, O, Sec for U is not part of the
       * Starndard matrix!
       */
      int numsymbol;
      /** 
       * minimum and maximum score in the matrix
       * This is a cache value to prevent repeated
       * calculation from the matrix.
       * maxs start with -1.
       * mins start with max_int
       */
      mutable int mins, maxs;
      /**
       * This field is populated only if allwords()
       * method is called.
       */
      mutable char **words;
      /**
       * capacity of words
       */
      mutable int wordsArraySize;
      /**
       * actual size of wordSize
       */
      mutable int wordSize;  // recoreds the size of words

      /**
       * Default installed matrix location.
       * should be /usr/local/share/orpara/matrix
       */
      static string default_path;

      /**
       * names of protein matrices in lower case
       */
      static const char* proteinMatrices[];
      /**
       * names of nucleic acid matrices in lower case
       */
      static const char* nucleicMatrices[];

   private:
      /**
       * Set both the mins(minimum) and maxs(maximum) from scanning the scoring
       * matrix. The mins and maxs is used for caching purpose only so that
       * repeated getter methods calls will not pose a performance problem.
       * This method is used by getMaxScore and getMinScore and should not be
       * used by any other method.
       */
      void setMinMaxScore() const;
};

/**
 * Scoring method for protein sequence comparison.
 * This class use aachar2num for look up. The code generated by Protein class
 * use the same function to ensure compatibility.
 */
class ProteinScoreMethod : public MatrixScoreMethod {
   public:
      /**
       * Default constructor use the default matrix and alphabet.
       */
      ProteinScoreMethod() 
         : MatrixScoreMethod() { 
            copyFrom(default_name, default_symb, 32, default_mat); }
      // the following is wrong because hashing of symbols
            //copyFrom(default_name, default_symb, default_numsymbol, default_mat); }

      ProteinScoreMethod(int go, int ge) 
         : MatrixScoreMethod(go, ge) {
            copyFrom(default_name, default_symb, 32, default_mat); }
            //copyFrom(default_name, default_symb, default_numsymbol, default_mat); }

      ProteinScoreMethod(const string& matrixName);

      ProteinScoreMethod(const string& matrixName, int go, int ge);
      /**
       * Copy constructor.
       */
      ProteinScoreMethod(const ProteinScoreMethod &psm) : MatrixScoreMethod(psm) { }
      ProteinScoreMethod& operator=(const ProteinScoreMethod &psm) {
         MatrixScoreMethod::operator=(psm); return *this; }

      /**
       * Overwrite base class method.
       */
      bool read(const string &p, const string &n);

      /**
       * Use a particular protein scoring matrix by name, such as blosum62.
       * @param matName is the name of the scoring matrix.
       * This is the same as the constructor(string).
       */
      void use(const string &matName);

      /**
       * use the aachar2num function in bioseq to
       * convert letters to numbers with A as zero, B as 1, ...
       * All letters must be upper case
       * Depends on the type of matrix, if nucleic acid matrix, then
       * use hashbase function,
       * else use aachar2num function.
       * This method will be implemented by derived class
       * For performance reason, we canot put it here.
       * the lookup function is used in tight loop
       * even the if/else statement is not allowed!
       */
      int lookup(char row, char col) const { 
         return mat[aachar2num(row)][aachar2num(col)]; 
      }
      /**
       * Overload parent class lookup method using
       * integer as input. If this method is not defined
       * then the char version will be used.
       */
      int lookup(int row, int col) const { 
         return MatrixScoreMethod::lookup(row, col);
      }

      void show(ostream &ous) const { MatrixScoreMethod::show(ous, &aachar2num); }
      void show() const { show(cout); }


      // static members
      /** blosum50 */
      static string default_name;
      static char default_symb[32];
      /** refer to the elements excluding the terminator in 
       * default_symb.  24 without terminator '\0' */
      static int default_numsymbol;
      /**
       * The default protein matrix.
       * This is blosum50. At construction time, this is used to initialize the
       * mat member from the parent class.
       */
      static int default_mat[32][32];
};

/**
 * Nucleic acid scoring method. 
 * hashbase function is used here.
 * The default matrix is NUC.4.4 without specifying any name.
 * You can also specify a particular matrix.
 * Under multi-threading environment, this class cause segmentation 
 * fault, but the ScoreMethod is fine.
 */
class NucleicScoreMethod : public MatrixScoreMethod {
   public:
      // the default_numsymbol needs to be checked!
      NucleicScoreMethod() 
         : MatrixScoreMethod() 
      { 
         copyFrom(default_name, default_symb, default_numsymbol, default_mat); 
      }

      NucleicScoreMethod(int go, int ge) 
         : MatrixScoreMethod(go, ge) 
      {
         copyFrom(default_name, default_symb, default_numsymbol, default_mat); 
      }

      /**
       * If a full path is given then it will make a new
       * matrix from it.
       */
      NucleicScoreMethod(const string& matrixName);
      /**
       * path to matrix and matrix name.
       * @param dir path to matrix directory. It is usualy
       *   installed in /usr/local/share/orpara/matrix.
       */
      NucleicScoreMethod(const string& dir, const string& matrixName);

      NucleicScoreMethod(const string& matrixName, int go, int ge);

      /**
       * Copy constructor uses the base class version.
       * This class does nothing.
       */
      NucleicScoreMethod(const NucleicScoreMethod &nsm) 
         : MatrixScoreMethod(nsm) 
      { 
//#ifdef DEBUG
         //cerr << "Called NucleicScoreMethod copy constructor\n";
         //show(cerr);
//#endif
      }
      NucleicScoreMethod(NucleicScoreMethod &&nsm) : MatrixScoreMethod(std::move(nsm)) { }
      /**
       * Assignmenet operator.
       */
      NucleicScoreMethod& operator=(const NucleicScoreMethod &nsm) {
         MatrixScoreMethod::operator=(nsm); return *this; }
      NucleicScoreMethod& operator=(NucleicScoreMethod &&nsm) {
         MatrixScoreMethod::operator=(std::move(nsm)); return *this; }

      /**
       * Read a nucleic acid matrix given path  p and name n
       * @param p a path to nucleic matrix directory
       * @param n matrix name.
       */
      bool read(const string &p, const string &n);

      /**
       * to use a particular Nucleic Acid scoring method by name.
       * @param method the string name of the method such as Nuc4.4
       * This function will set the internal matrix to be used for scoring.
       * The scoring matrix must be stored under /path member.
       */
      void use(const string &matName);

      /**
       * Explicity call parent class method to prevent char to int implicity
       * conversion by the compiler.
       */
      int lookup(int row, int col) const { 
         return MatrixScoreMethod::lookup(row, col);
      }

      /**
       * Use the hashbase function for nucleic acids.
       */
      int lookup(char row, char col) const { 
         return mat[hashbase(row)][hashbase(col)]; 
      }

      void show(ostream &ous) const { MatrixScoreMethod::show(ous, &hashbase); }
      void show() const { show(cout); }

      ////// static members //////////
      /** NUC.4.4 */
      static string default_name;
      static char default_symb[32];
      /** 24 with terminator '\0' */
      static int default_numsymbol;
      /**
       * The default nucleic acid matrix.
       * This is the NUC.4.4 
       */
      static int default_mat[32][32];
};
}

#endif
