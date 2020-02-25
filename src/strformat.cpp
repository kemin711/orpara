#include "strformat.h"
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cassert>
#include <string>
#include <stdexcept>
#include <cmath>

//#include <algorithm>
//file: strformat.cpp
//
using namespace std;
//using namespace KZUtility;

namespace orpara {
/**
 * efficient for vectors with less than a few hudrend elements
 * large sets needs to use the set algorithm
 * @return true of two containers share common elements.
 */
bool intersect(vector<string> &s1, vector<string> &s2) {
	for (unsigned int i=0; i<s1.size(); i++) {
		for (unsigned int j=0; j < s2.size(); j++) {
			if (s1[i] == s2[j]) return true;
		}
	}
	return false;
}

/**
 * @return the number of words in the string
 */
int wc(const string &str) {
	string::size_type i = str.find_first_not_of(' ');
	int count = 0;
	while (i != string::npos) {
		count++;
		i = str.find_first_of(' ', i+1);
		if (i == string::npos) break;
		i = str.find_first_not_of(' ', i+1);
	}
	return count;
}

bool isnumber(const string &str) {
   if (str.empty()) {
      cerr << "empty string is not a number!\n";
      return false;
   }
   string::size_type i = 0;
	while (isdigit(str[i])) i++;
	return i == str.size();
}

bool isnumber(const string &str, const string dc) {
   if (str.empty()) {
      return false;
   }
   string::size_type i = 0;
	while (i<str.size()) {
		if (!isdigit(str[i]) && dc.find(str[i]) == string::npos) 
			return false;
		i++;
	}
	return true;
}

bool isInt(const char str[]) {
   size_t i=0;
   while (str[i] != '\0') {
      if (!isdigit(str[i])) return false;
      ++i;
   }
   return true;
}

bool isupper(const string &str) {
   string::size_type i = 0;
	while (i<str.length()) {
		if (isalpha(str[i]) && !std::isupper(str[i]))
			return false;
		++i;
	}
	return true;
}

bool startwith(const string &str, const string &head) {
   return str.substr(0, head.length()) == head;
}

bool endwith(const string &str, const string &tail) {
   string::size_type i = str.rfind(tail);
   if (i != string::npos) {
      if (i+tail.length() == str.length()) 
         return true;
   }
   return false;
}

void upper(const char lo[], char *up)
{  // this is a helper function that convers string lo into all uppercase 
	// returns upper case
	int i = 0;
	while (lo[i] != '\0')  {
		up[i] = toupper(lo[i]);
		i++;
	}
	up[i] = '\0';
}

/* return a new string with lowercases */
string getLower(const string &str) {
	string tmp = str;
	for (string::size_type i=0; i<tmp.size(); i++) 
		tmp[i] = tolower(tmp[i]);
	return tmp;
}
void strTolower(string &str) {
	for (string::size_type i=0; i<str.size(); i++) 
		str[i] = tolower(str[i]);
}
void strTolower(const string &str1, string &str2) {
	str2=str1;
	for (string::size_type i=0; i<str2.size(); i++) 
		str2[i] = tolower(str2[i]);
}

void strToupper(string &str) {
	for (string::size_type i=0; i<str.length(); i++)
		str[i] = toupper(str[i]);
}

string lc(const string& str) {
	string tmp(str);
	string::iterator i = tmp.begin();
	while (i != tmp.end()) {
		*i = tolower(*i);
		++i;
	}
	return tmp;
} 

string uc(const string& str) {
   string tmp(str);
   for (string::size_type i=0; i<tmp.size(); ++i) {
      if (!std::isupper(tmp[i]))
         tmp[i]=toupper(tmp[i]);
   }
   return tmp;
}

int cmp_nocase(const string &s1, const string &s2) {
	string::const_iterator i1 = s1.begin();
	string::const_iterator i2 = s2.begin();

	while (i1 != s1.end() && i2 != s2.end()) {
		if (tolower(*i1) != tolower(*i2)) 
			return (tolower(*i1)<tolower(*i2))? -1 : 1;
		++i1;
		++i2;
	}
	if (s1.size() == s2.size()) return 0;
	else return (s1.size()<s2.size())? -1 : 1;
}

void rmsp1(const char a[], char b[])
{
	char *x = new char[strlen(a) + 1];
	strcpy(x, a);
	char *aa = x;
	int j = 0;
	while (*aa != '\0') {  //converts all \n to space
		if (*aa == '\n') *aa = ' ';
		aa++;
	}
	aa--;
	while ( isspace(*aa) ) aa--;
	*(aa+1) = '\0'; //this removes extra space from the end
	//the following code removes space from the middle
	aa = x;
	while ( *aa != '\0' )  {
		while ( !isspace(*aa) && *aa != '\0'  ) {
			b[j] = *aa;
			aa++;
			j++;
		}
		if ( *aa != '\0' ) {
			b[j] = ' ';
			aa++;
			j++;
			while ( isspace(*aa) ) aa++;
		}			
	}
	b[j] = '\0';
	delete[] x;
}

void rmsp(char a[])
{ //remove spaces from the middle and the ends
	if (!strstr(a, " ") && !strstr(a, "\n")) return;
	rmtsp(a); //remove trailing space
	
	//the following code removes space from the middle
	char *aa = a, *bb;
	while (*aa != '\0') {  //converts all \n to space
		if (*aa == '\n') *aa = ' ';
		aa++;
	}
	
	aa = a;
	while ( !isspace(*aa) && *aa != '\0'  ) {
		aa++;
	}
	if ( *aa == '\0' ) return;
	if ( isspace(*aa) ) { //not the end yet, found the space
		aa++;
		bb = aa;
		while ( isspace(*aa) ) aa++;
	}
	while ( *aa != '\0' )  {
		while ( !isspace(*aa) && *aa != '\0'  ) {
			*bb = *aa;
			aa++;
			bb++;
		}
		if ( *aa != '\0' ) {
			aa++;
			*bb = ' ';
			bb++;
			while ( isspace(*aa) ) aa++;
		}
	}
	*bb = '\0';
}

void rmsp(string& str) {
   str.erase(remove_if(str.begin(), str.end(), [](char c){ return isspace(c); }), str.end());
}

// not tested yet
void removeWhiteChar(string& str) {
   string::size_type i=0,j=0;
   while (j<str.size()) {
      if (i<j) { str[i]=str[j]; }
      if (isspace(str[j])) {
         ++j;
      }
      else {
         ++i; ++j;
      }
   }
   if (i < str.size()-1) str.resize(i);
}

/**
 * @return a new string with all c deleted in str
 */
string delall(const string& str, const char c) {
	string tmp;
   tmp.reserve(str.size());
   for (const char x : str) {
      if (x != c) tmp += x;
   }
	return tmp;
}

/**
 * Same as delall
 */
string deleteChr(const string &str, char c) {
   string tmp;
   for (string::size_type i=0; i<str.size(); i++) {
      if (str[i] != c) {
         tmp.push_back(str[i]);
      }
   }
   return tmp;
}


string delall(const string& str, const string& c) {
	string tmp = str;
   string::size_type i;
	i = tmp.find_first_of(c);
	while (i != string::npos) {
		tmp.erase(i, 1);
		i = tmp.find_first_of(c, i);
	}
	return tmp;
}

void rmdquote(char a[])
{
   char *xx = strchr(a, '\"');
   while (xx) {  //replacing " with space!!
      *xx = ' ';
      xx = strchr(++xx, '\"');
   }
}
void rpldquote(char a[], const char r)
{
	char *x = strchr(a, '\"');
	while (x) {
		*x = r;
		x = strchr(++x, '\"');
	}
}

void removeQuote(string& str, const char q) {
   if (str[0] == q && str.back() == q) {
      str = str.substr(1, str.size()-2);
   }
}

string getQuoteLess(const string& str, const char q) {
   if (str[0] == q && str.back() == q) {
      return str.substr(1, str.size()-2);
   }
   if (str[0] == q || str.back() == q) {
      throw runtime_error(string(__func__) + ":ERROR: unbalanced quote: " + str);
   }
   return str;
}

void removeDoubleQuote(string& str) {
   if (str[0] == '"' && str.back() == '"') {
      str = str.substr(1, str.size()-2);
   }
}

string tr(const string& str, char i, char o) {
	string::const_iterator it = str.begin();
	string tmp;
	tmp.reserve(str.size());
	while (it != str.end()) {
		if (*it == i) tmp += o;
		else tmp += *it;
		it++;
	}
	return tmp;
}

bool isAlpha(const string& str) {
   for (size_t i=0; i<str.size(); ++i) {
      if (!isalpha(str[i])) return false;
   }
   return true;
}

int isName(char *n)
{
	if (strstr(n, " ")) {
		if (strstr(n, "virus") || strstr(n, "Virus")) 
			return 1;
		else {
			if (std::isupper(n[0])) {
				char *p = n + 1;
				int numwd = 0;
				while (*p != '\0' && numwd < 2) {
					if (std::isupper(*p)) return 0;
					if (isspace(*p)) numwd++;
						p++;
				}
				return 1;
			}
			else return 0;
		}
	}
	else return 0;
}

void singleSpace(string &str) {
	string::size_type i = 0;
	string::size_type s;

	while (i<str.size()) {
		if (isspace(str[i])) {
			s = i;
			++i;
			while (i<str.size() && isspace(str[i])) {
				i++;
			}
			if (i>s+1) str.erase(s, i-s-1);
		}
		else {
			while (i<str.size() && !isspace(str[i])) i++;
		}
	}
}
void rmtsp(char a[])
{
	char *p = a + strlen(a) - 1;
	while ( isspace(*p) || *p == '\n' ) p--;
	*(p+1) = '\0';
}
	
void newline(istream &ins)
{
	char ch;
	ins.get(ch);
	while (ch != '\n' && ch != '\r') ins.get(ch);
}

void dlc(char a[])
{
	if (strlen(a) > 0) {
		a[strlen(a) -1] = '\0';
	}
	else cout << "string has no character in it\n";
}

int getNumber(char *&ptr)
{ //get number from the location format, advance the char pointer
	//to the next non-digit character.  The pointer must be set to
	//the first digit			
	char number[11];
	 int i = 0;
	 while (isdigit(*ptr)) {
		 number[i] = *ptr;
		 i++;
		 ptr++;
	 }
	 number[i] = '\0';
	 return (atoi(number));
} 
				
/**
 * Collect all integer numbers from str
 */
vector<int> getAllInt(const string &str) {
   string::size_type i=0;
   string::size_type b;
   vector<int> tmp;
   while (i<str.size()) {
      while (i<str.size() && !isdigit(str[i])) ++i;
      if (i >= str.size()) break;
      b=i;
      while (i<str.size() && isdigit(str[i])) ++i;
      tmp.push_back(atoi(str.substr(b,i-b).c_str()));
   }
   return tmp;
}

int getInt(const string &str) {
   string::size_type i=0;
   while (i<str.size() && !isdigit(str[i])) ++i;
   if (i==str.size()) {
      cerr << "string " << str << " has no digit!";
      exit(1);
   }
   string::size_type j=i+1;
   while (j<str.size() && isdigit(str[j])) ++j;
   return atoi(str.substr(i,j-i).c_str());
}
				
int itoa(unsigned int n, char a[])
{  //returns the number of digits.  the size of string A must be large enough
	// 10 is good enough for most numbers
	//32-bit unsigned = 2x 2147483648.  10 digits.    
	if (n == 0) {
		a[0] = '0';
		a[1] = '\0';
		return 1;
	} //convets 0 digit to 0 string
	//I used to convets 0 digit to null string

	char s[20]; //20 give you max on all computers
	int i = 0, j=0;
	while (n > 0)  {
		s[i++] = (char)( (int)'0' + n%10 );
		n = n/10;
	}
	s[i] = '\0';
	//s string is in reversed order.  We need to reverse it
	i--;
	while (i>=0) {
		a[j]=s[i];	
		j++;
		i--;
	}
	a[j]='\0';
	return j;
}

/* efficiently translate integer into string 
 * using the stringstream object
 * */
string itos(int n) {
	ostringstream ostr;
	ostr << n;
	return ostr.str();
}

string ftos(float f, int nd) {
	ostringstream ostr;
   int tens = stoi("1" + string(nd, '0'));
	ostr << round(f*tens)/tens;
	return ostr.str();
}

/* here is two other ways of converting integer to string
both use recursive functions
void itoa(int n, char a[])
{
	static int i=0;
	if (n<10) {
		cout << char(n);
		a[i]=char( int('0') + n);
		i++;
		a[i] = '\0';
	}
	else {
		a[i] = char( int('0') + n%10);
		i++;
		n = n/10;
		itoa(n, a);
	}
}

void itoa(int n)
{
	static int i=0;
	if (n<10) {
		cout << char( int('0') + n);
	}
	else {
		cout << char ( int('0') + n%10);
		n = n/10;
		itoa(n);
	}
}
*/

int substr(const char ln[], int s, char sub[])
{ //s is the start index, base 0
	const char *p = ln + s;
	if (isspace(*p) || *p == '\0') {
		sub[0]='\0';
		return 0;
	}
	else {
		const char *pp = p+1;
		while (!isspace(*pp) && *pp != '\0') pp++;
		int i = pp-p;
		if (*pp == '\0')  strcpy(sub, p);
		else {
			strncpy(sub, p, i);
			sub[i]='\0';
		}
		return i;
	}
}

int substr(const char ln[], int s, int f, char sub[])
{ //s is the start index, f is the finishing index, base 0
	const char *p = ln + s;
	const char *pp = ln + f;
	while (isspace(*p) && p < pp)  p++;  //skipping whitespace
	if (p == pp) { // empty sub string string
		sub[0] = '\0';
		return 0;
	}

	while (isspace(*pp)) pp--;
	int i = pp-p+1;
	memcpy(sub, p, i);
	sub[i]='\0';
	return i;
}

void firstwd(const char ln[], char sub[], char term)
{
	int i = 0;
	while (ln[i] != term && ln[i] != '\0') i++;
	strncpy(sub, ln, i);
	sub[i] = '\0';
}
string firstword(const string &str, const string &delim) {
	string::size_type i = str.find_first_not_of(delim);
	string::size_type b = str.find_first_of(delim, i);
	return str.substr(i, b-i);
}
void lastwd(const char ln[], char wd[], char sep) {
	int i = strlen(ln) - 1;
	while (i>0 && ln[i] != sep) i--;
	if (i != 0) i++;
	strcpy(wd, ln+i);
}
string lastword(const string &str, const string &delim) {
	string::size_type i = str.find_last_not_of(delim);
	string::size_type b = str.find_last_of(delim);
	return str.substr(b+1, b-i);
}
string acronym(const string &str, int n) {
	string tmp;
	string::size_type i=0;
	int l;
	while (i<str.size() && !isalpha(str[i])) i++; // first word char
	while (i < str.size()) {
		l=0;
		while (l<n && i<str.size()) {
			if (isalpha(str[i])) {
				tmp += str[i++];
				l++;
			}
			else if (str[i] == '-') ++i;
			else break;
		}
		while (i<str.size() && (isalpha(str[i]) || str[i] == '-')) {
			// skip the rest of the word
			i++;
		}
		while (i<str.size() && !isalpha(str[i])) i++; 
	}
	return tmp;
}
string acronymWithDigit(const string &str, int n) {
	string tmp;
	string::size_type i=0;
	int l;
	while (i<str.size() && !isalpha(str[i]) && !isdigit(str[i])) 
		i++; // first word char
	while (i < str.size()) {
		l=0;
		while (l<n && i<str.size()) {
			if (isalpha(str[i]) || isdigit(str[i])) {
				tmp += str[i++];
				l++;
			}
			else break;
		}
		while (i<str.size() && (isalpha(str[i]) || isdigit(str[i]))) {
			i++;  // rest of the word
		}
		while (i<str.size() && !isalpha(str[i]) && !isdigit(str[i])) 
			i++; 
	}
	return tmp;
}

string acronymWithTag(const string&str, int n) {
	vector<string> tmpv = dissect(str, " ,");
   string::size_type i;
	string tmp;
	for (i=0; i<tmpv.size(); i++) {
		if (tmpv[i].length() > (unsigned int)n) tmp += tmpv[i].substr(0,n);
		else tmp += tmpv[i];
	}
	string lwd = tmpv[tmpv.size()-1];
	i = lwd.length()-1;
	if ((unsigned int)n < lwd.length() && isdigit(lwd[i])) {
		--i;
		while (i>0 && (isdigit(lwd[i]) || lwd[i] == '-'
				|| lwd[i] == '.')) --i;
		++i; // pointing to the digit or -
		if (i >= (unsigned int)n) tmp += lwd.substr(i);
	}
	return tmp;
}

void append(char *&head, const char *tail, int &len, int &maxlen, 
		int incr)
{
	int tl = strlen(tail);
	if (tl > incr) incr = tl;
	char *temp;
	
	if ((len + tl) > maxlen) {
		maxlen += incr;
		temp = new char[maxlen+1];
		//assert(temp != 0);
		memcpy(temp, head, len);
		temp[len]='\0';
		delete[] head;
		head = temp;
	}
	strcat(head, tail);
	len += tl;
}

void append(char *&head, const char *tail, int &len, int &maxlen)
{
	int tl = strlen(tail);
	char *temp;
	if ((len + tl) > maxlen) {
		maxlen *= 2;
		temp = new char[maxlen+1];
		//assert(temp != 0);
		memcpy(temp, head, len);
		temp[len]='\0';
		delete[] head;
		head = temp;
	}
	strcat(head, tail);
	len += tl;
}

/****** helper function for wjoinseg ************/
vector<string> split(const string &str, const char sep) {
	vector<string> tmp;
	string::size_type i=0,ii;
	while ((ii=str.find(sep, i)) != string::npos) {
		tmp.push_back(str.substr(i,ii-i));
		i= ii+1;
	}
	if (i == str.length()) tmp.push_back("");
	else tmp.push_back(str.substr(i));
	return tmp;
}

// optional quote
// will preserve the quote in the value
vector<string> splitQuoted(const string &str, const char quote, const char sep) {
	string::size_type i=0,ii, iq;
   bool quoted = false;
   if (isspace(str[0])) {
      cerr << "WARN: string starts with white space!\n";
      ++i;
   }
   if (str[i] == quote) {
      quoted=true;
   }
	vector<string> tmp;
   do {
      if (quoted) { // use quote
         iq = str.find(quote, i+1);
         if (iq == string::npos) {
            cerr << "filed missing right quote!\n";
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR: field missing right quote: " + str);
         }
         tmp.push_back(str.substr(i, iq-i+1));
         quoted = false;
         i = iq + 1; // either quote or end of line
         if (i < str.size()) {
            assert(str[i] == sep);
            ++i;
         }
      }
      else { // not quoted string field<SEP>abcdef<SEP>
         //5,7,,,8,ABC => 5|7|||8|ABC
         if (i < str.size() && str[i] == sep) {
            tmp.push_back(string());
            ++i;
         }
         /*
         if (i+1 < str.size() && str[i+1] == sep) { // special empty field case
            tmp.push_back(string());
            ++i;
         }
         */
         else {
            ii = str.find(sep, i);
            if (ii == string::npos) {
               tmp.push_back(str.substr(i));
               break;
               //i=ii;
            }
            else {
               tmp.push_back(str.substr(i,ii-i));
               i = ii + 1;
            }
         }
      }
      if (i != str.size() && i != string::npos && str[i] == quote) {
         quoted = true;
      }
	}
   while (i != str.size() && i != string::npos);
   if (str.back() == sep) { // last field is empty, must append empty string!
      tmp.push_back(string());
   }
	return tmp;
}

/* sep is a separator such as ',' or "..".  the location string can be used
 *  * as input
 *   * */
vector<string> split(const string &str, const char sep[]) {
	string::size_type i=0,j;
	int seplen = strlen(sep);  // length of separator
	vector<string> vec;
	while ((j=str.find(sep, i)) != string::npos) {
		vec.push_back(str.substr(i, j-i));
		//cout << str.substr(i, j-i) << " testing\n";
		i = j+seplen;
	}
	vec.push_back(str.substr(i));
	return vec;
}
// use any one character in delim to split the string
vector<string> splitOneOf(const string &str, const string &delims) {
	string::size_type begIdx=0, endIdx;
	vector<string> words;

	endIdx = str.find_first_of(delims);
	//begIdx = str.find_first_not_of(delims);
	while (endIdx != string::npos) {
		words.push_back(str.substr(begIdx, endIdx-begIdx));
		begIdx = endIdx + 1;
		endIdx = str.find_first_of(delims, begIdx);
	}
	if (begIdx>str.size()) words.push_back("");
	else words.push_back(str.substr(begIdx));
	return words;
}
vector<string> dissect(const string &str, const string &delims)
{
	string::size_type begIdx, endIdx;
	vector<string> words;

	begIdx = str.find_first_not_of(delims);
	while (begIdx != string::npos) {
		endIdx = str.find_first_of(delims, begIdx);
		if (endIdx == string::npos) {
			words.push_back(str.substr(begIdx));
			return words;
		}
		else words.push_back(str.substr(begIdx, endIdx-begIdx));
		begIdx = str.find_first_not_of(delims, endIdx);
	}
	return words;
}
set<string> digest2set(const string& str, const string& delims) {
	string::size_type i, ii;
	set<string> words;
	i = str.find_first_not_of(delims);
	while (i != string::npos) {
		ii = str.find_first_of(delims, i);
		if (ii == string::npos) {
			words.insert(str.substr(i));
			return words;
		}
		else words.insert(str.substr(i, ii-i));
		i = str.find_first_not_of(delims, ii);
	}
	return words;
}

pair<string,string> breakString(const string &str, const string &sep) {
   string::size_type i = str.find(sep);
   if (i == string::npos) {
      cerr << "string " << str << " has no separator: '"
            << sep << "' returning a pair of empty strings.\n";
      return pair<string,string>("", "");
      //exit(1);
   }
   return make_pair(str.substr(0,i), str.substr(i + sep.length()));
}

string join(const vector<string>& parts, char sep) {
   string tmp = parts[0];
   for (size_t i=1; i<parts.size(); ++i) {
      tmp += (sep + parts[i]);
   }
   return tmp;
}

void trim(string &str) {
	string::size_type idx = str.find_last_not_of(" \t\n\v\f\r");
	str.erase(idx+1);
}

/**
 * remove space in front of str
 */
void trimLeadingSpace(string &str) {
   string::size_type i=0;
   while (i<str.size() && isspace(str[i])) ++i;
   str=str.substr(i);
}

string trimSpace(const string &str) {
   if (str.empty()) return str;
   string::size_type i,j; 
   i=0;
   while (i<str.size() && isspace(str[i])) ++i;
   if (i == str.size()) 
      return string(); // empty string
   j = str.size() - 1;
   //while (j>0 && isspace(str[j])) --j;
   while (isspace(str[j])) --j; // above code ensures j>0
   return str.substr(i, j-i+1);
}

string chopend(const string &str, const string &sub) {
   if (sub.length() >= str.length()) 
      return str;
   if (str.substr(str.length()-sub.length()) == sub) {
      return str.substr(0, str.length()-sub.length());
   }
   else {
      return str;
   }
}

string chopfront(const string &str, const string &sub) {
   if (sub.length() >= str.length()) 
      return str;
   if (str.substr(0, sub.length()) == sub) {
      return str.substr(sub.length());
   }
   else return str;
}

string chopLastWord(const string &str) {
   string::size_type i = str.rfind(' ');
   if (i != string::npos) {
      return str.substr(0,i);
   }
   return str;
}

string chopFirstWord(const string &str) {
   string::size_type i = str.find(' ');
   if (i != string::npos) {
      return str.substr(i);
   }
   return str;
}

void writeSequence(const string &seq, ostream &ous, const int width) 
{
   string::size_type i=0;
	while (i<seq.length()) {
		ous << seq.substr(i, width) << endl;
		i += width;
	}
}

// using a stringstream would be simplier to implement
vector<int> extractInt(const string &str) {
	string::size_type i=0, j;
	vector<int> tmp;

	while (i<str.size() && !isdigit(str[i])) i++;
	while (i<str.size()) {
		j=i+1;
		while (j<str.size() && isdigit(str[j])) j++;
		tmp.push_back( atoi(str.substr(i,j-i).c_str()) );
		i=j+1;
		while (i<str.size() && !isdigit(str[i])) i++;
	}
	return tmp;
}

string str2upper(const string &str) {
   string tmp(str);
	for (string::size_type i=0; i<tmp.length(); i++)
		tmp[i] = toupper(tmp[i]);
   return tmp;
}

bool codonIsStop(const string &codon) { 
   string tmp(codon);
   strToupper(tmp);
   return tmp == "TAA" || tmp == "TAG" || tmp == "TGA"; 
}

bool subseqIsStop(const string &seq, int i) { 
   string tmp(seq.substr(i,3));
   strToupper(tmp);
   return tmp == "TAA" || tmp == "TAG" || tmp == "TGA"; 
}

bool isScientificSpeciesName(const string &str) {
   if (str.find("uncultured") != string::npos ||
         str.find("unidentified") != string::npos)
      return false;
   vector<string> words = split(str, ' ');
   if (words.size() < 2) return false;
   if (words[1] == "bacterium" || words[1] == "sp."
         || words[0] == "human") return false;
   if (startwith(str, "Incertae sedis")) {
      if (words.size() < 4) {
         return false;
      }
      else {
         words.erase(words.begin(), words.begin()+2);
      }
   }
   if (std::isupper(words[0][0]) && islower(words[1][0])
      && isAlpha(words[0]) && isAlpha(words[1])) {
      return true;
   }
   else if (words[0] == "Candidatus" && words.size() == 3
         && std::isupper(words[1][0]) && islower(words[2][0])) {
      return true;
   }
   return false;
}

string getScientificSpeciesName(const string &taxon) {
   vector<string> words=split(taxon, ' ');
   if (std::isupper(words[0][0]) && islower(words[1][0])
         && isAlpha(words[0]) && isAlpha(words[1]))
   {
      return words[0] + " " + words[1];
   }
   else if (words.size() > 2) {
      if (words[0] == "Candidatus" && std::isupper(words[1][0])) {
         return words[0] + " " + words[1] + words[2];
      }
      else {
         throw runtime_error("ERROR: " + taxon + " does not contain species name");
      }
   }
   else {
      cerr << "ERROR: " << taxon << " is not a proper species name\n";
      exit(1);
   }
   return taxon;
}

string fileBasename(const string& pathstr) {
   string::size_type x = pathstr.rfind('/');
   string result = pathstr;
   if (x != string::npos) {
      result = pathstr.substr(x+1);
   }
   return result;
}

string fileBasename(const string& pathstr, const string& suffix) {
   string base = fileBasename(pathstr);
   string::size_type x = base.rfind(suffix);
   if (x != string::npos) {
      base = base.substr(0,x);
   }
   return base;
}

string getFileStem(const string& filename) {
   string tmp = fileBasename(filename);
   string::size_type i = tmp.rfind('.');
   if (i != string::npos) {
      return tmp.substr(0,i);
   }
   return tmp;
}


// end of orpara namespace
}
