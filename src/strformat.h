// (c) 1997 Kemin Zhou at The Molecular Sciences Institute
// strformat.h
#ifndef STRFORMAT_H
#define STRFORMAT_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <set>
#include <cstring>

//this file is upcase.h 
//namespace KZUtility {
// causing too much trouble

using namespace std;
//////// Information about string //////////

namespace orpara {
/**
 * count the number of words
 */
int wc(const string &str);  // word cound
bool isnumber(const string &str); // all digits
// whether is it a number with other chars
// such as 1.2, 1,000, -1
bool isnumber(const string &str, const string dc);
/** 
 * @return true if all of the chars are in upper case
 */
bool isupper(const string &str);
/**
 * All chars are alpha
 */
bool isAlpha(const string& str);
bool endwith(const string& str, const string &tail);
bool startwith(const string& str, const string &head);

/** 
 * Case transformation function 
 */
void upper(const char lo[], char *up);
/**
 * get an all lower case version of the original string.
 */
string getLower(const string &str);
void strTolower(string &str);
void strTolower(const string &str1, string &str2);
/**
 * @param str this string will become upper case after the function call.
 */
void strToupper(string &str);
string str2upper(const string &str);
string lc(const string &str);

/** compare to string and ignoring case **/
int cmp_nocase(const string &s1, const string &s2);
bool intersect(vector<string> &s1, vector<string> &s2);

/****** Modification  *****
 * replace, deletion
 */
/** delete all character from a string
 * @see also deleteChr
 * */
string delall(const string& str, const char c);
/** delete all string c from str */
string delall(const string& str, const string& c);
/** alias for delall
 * @return a string with all 'c' removed from it. 
 * if 'c' not present in the input, the 
 * return the original string.
 * */
string deleteChr(const string &str, char c);

/** converting all character i to o in str */
string tr(const string& str, char i, char o);
/** removing leading and trailing spaces or other while char if exist */
void trim(string &str); 
/**
 * remove non-printable char from both ends
 */
string trimSpace(const string &str); 
/** remove all leading spaces */
void trimLeadingSpace(string &str);

/**
 * chopend removes substring from tail
 * if str end with sub, then remove the sub and return
 * the shorter string.
 * If str does not end with sub, then return the 
 * original string.
 */
string chopend(const string &str, const string &sub);
string chopfront(const string &str, const string &sub);
string chopLastWord(const string &str);
string chopFirstWord(const string &str);

/** delete last character from an C-string a */
void dlc(char a[]); //deltes the last character
void rmtsp(char a[]); //remove trailing white characters such as space and \n
void rmsp(char a[]);
//removes extra white-char (space, \t, and \n) from the string, so 
//that only one space is left to separate words
void singleSpace(string &str);
/** remove all spaces */
void rmsp(const char a[], char b[]);
//overloaded version.   input string a, output string b
void rmdquote(char a[]); //removes double quote
void rpldquote(char a[], const char r); //replace double quote with r

/** This function is not so good and should be removed in the future
 */
void newline(istream &ins); //removes remaining character from input stream

//////////////////////////////////////////////////////////

/** Conversion to other types such as integer, double
 */
int getNumber(char *&ptr);  
/** return the first integer number in the string */
int getInt(const string &str);
/** extract all integer numbers in the string 
 * @param str string contain multiple integers
 *    that can be separate by any number of any non-digit
 *    characters.
 * */
vector<int> getAllInt(const string &str);


//get number from a pointer that is on the first digit of a string.  
//The pointer will be advanced to the first non-digit character
/** reverse of atoi
 * convert into an array of char
 */
int itoa(unsigned int n, char a[]);
/** convert integer type to string type
 * the toString() method is more generic
 */
string itos(int n);
/**
 * converting floating number to string.
 * @param f input floating number
 * @param nd number of digits after the point
 */
string ftos(float f, int nd);

/** this can be used to conver most of the build-in types
 * to string
 */
template<class T> string toString(const T &val) {
   ostringstream ous;
   ous << val;
   return ous.str();
}

/** Make a long name so it will be different
 * I tried to use name space but I have a hard time
 * to compile and link!
 */
template<class T> string anyToString(const T &val) {
   ostringstream ous;
   ous << val;
   return ous.str();
}
vector<int> extractInt(const string &str);

////////////////////////////////////////////////////////////

/** Extraction */

/** picks substring from index s to a white space of ln, copy to sub
 * sub is set to the substring without any whitespace at beginning
 * or end. If from s to the end is all whitespace, then return 0
 * and sub is set to empty C string.
 * returns strlen of the substr, 0 if no substr, or substr is empty
 * s is 0-based index
 */
int substr(const char ln[], int s, char sub[]);

/** picks non-whitespace substring from s to f
 * returns the length of the substring
 * s starting index, f ending index (inclusive)
 */
int substr(const char ln[], int s, int f, char sub[]);

void firstwd(const char ln[], char sub[], char term=' ');
/**
 * extract the first word of the string
 */
string firstword(const string &str, const string &delim=" ,.()");
void lastwd(const char ln[], char wd[], char sep = ' ');
/**
 * @return the last word from a string str.
 */
string lastword(const string &str, const string &delim=" ,.()");

/** acronym based on, at most, n letters of each word
 */
string acronym(const string&str, int n=1);
string acronymWithDigit(const string&str, int n=1);
// if taild with digit, then all will be added
string acronymWithTag(const string&str, int n=1);

void append(char *&head, const char *tail, int &len, int &maxlen, 
		int incr);
//carefull, these are managed by new allocation
//incr is the size increment, should be larger than strlen(tail)
//if not the function will increase it automatically
void append(char *&head, const char *tail, int &len, int &maxlen);
//////////////////////////////////////////////////////////////////

/******** Digestion ***********/

/** A version of split that will separate a string into elements
 *  separated by sep that is a single character.
 *  It even returns empty string as elements if 
 *  nothing is in between sep. 
 *  @param sep is the separator as single char.
 */
vector<string> split(const string &str, const char sep='\t');

/** separate str into an array of strings using sep as delimiter sep is a
 * string Same behavior as the char version seprator.  Anything in betwee sep
 * will be an array element.  For example, if sep is ??  then ??tabc??  will
 * have three elements, the first and the third are empty strings
 *
 * @param sep is the delimiter that is a c_string as compare
 *    to the single character version.
 */
vector<string> split(const string &str, const char sep[]);

/** use one of the separators ( \t,.;), no more
 * will discard the leading separators
 * This function should be removed, it is confusing
 * All split functions should not discard leading fields.
 * Needs to redefine.
 */
//vector<string> split(const string &str); //default separator ( \t,.;)

// will split the string use any one char in delims
vector<string> splitOneOf(const string &str, const string &delims = ",. \t");

/** anything in delims will be discarded
 * only the none-delims char are left, and packed in the vector
 * more or less digestions
 * same as split, but will use any of the characters as separators
 * leading and trailing separators are discarded
 */
vector<string> dissect(const string &str, const string &delims = ",. \t");
set<string> digest2set(const string &str, const string &delims=",. \t");

/** break a string into two parts according the the 
 *  separator sep
 *  @return a pair of strings.  if sep is not found inside str, 
 *    then I will return two empty strings as a pair.
 */
pair<string,string> breakString(const string &str, const string &sep);

//////////////////////////////////////////////////////////////////////////

//write sequence to ous line by line; seq is a long string
void writeSequence(const string &seq, ostream &ous, const int width=70);

/**
 * is scientific name
 */
int isName(char *n);

///////////////////////// Bioinformatics Helpers /////////////////
bool codonIsStop(const string &codon);
//{ return codon == "TAA" || codon == "TAG" || codon == "TGA"; }
//bool codonIsStop(char codon[3]) { return !strcmp(codon, "TAA") || !strcmp(codon, "TAG") || !strcmp(codon, "TGA"); }
/** seq must be in upper case
 * @param i. is the index of the subsequence of 3 nt.
 **/
bool subseqIsStop(const string &seq, int i);
//{ return seq.substr(i,3) == "TAA" || seq.substr(i,3) == "TAG" || seq.substr(i,3) == "TGA"; }
/**
 * To check that the string is a scientific name
 * For example Sacchromyces cerevisiae
 * Subspecies and strain lableling can follow the species.
 */
bool isScientificSpeciesName(const string &str);
/**
 * Extract the species name without any strain labeling.
 */
string getScientificSpeciesName(const string &taxon);
}

#endif
