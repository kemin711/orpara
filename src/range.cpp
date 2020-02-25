#include <sstream>
#include <vector>
#include <functional>
#include <limits>
#include <cstdlib>
#include <utility>

#include "range.h"
#include "strformat.h"

namespace orpara {

int Range::margin=16;
float Range::marginFraction=0.1;

Range& Range::operator=(const Range &r) {
	if (this != &r) {
		b=r.b;
		e=r.e;
	}
	return *this;
}

ostream& operator<<(ostream &ous, const Range &r) {
   ous << r.b << "->" << r.e;
   return ous;
}

// if merge successful return true
// else return false.
bool Range::merge(const Range &r, bool sd) {
//int Range::merge(const Range &r) {
	// only merge if overlap
	if (overlap(r,sd) == 0) return false;

   if (direction() == '+' && r.direction()== '+') {
      b=b<r.b? b : r.b;
      e=e>r.e? e : r.e;
   }
   else if (direction() == '-' && r.direction() == '-') {
      b=b>r.b? b : r.b;
      e=e<r.e? e : r.e;
   }
   else {
      if (!sd) {
         // how to merge fragments of different direction?
         // The resulting Range is always in the '+' direction
         int bb = smallerEnd() < r.smallerEnd() ? smallerEnd() : r.smallerEnd();
         int ee = largerEnd() > r.largerEnd() ? largerEnd() : r.largerEnd();
         b=bb;
         e=ee;
         return true;
      }
      else return false;
   }
   return true;
}

int Range::combine(const Range &r) {
   int olp=overlap(r);
   if (direction() == '+' && r.direction()== '+') {
      b=b<r.b? b : r.b;
      e=e>r.e? e : r.e;
   }
   else if (direction() == '-' && r.direction() == '-') {
      b=b>r.b? b : r.b;
      e=e<r.e? e : r.e;
   }
   else {
      // how to merge fragments of different direction?
      // The resulting Range is always in the '+' direction
      // operation modified object cause big trouble!
      int bb = smallerEnd() < r.smallerEnd() ? smallerEnd() : r.smallerEnd();
      int ee = largerEnd() > r.largerEnd() ? largerEnd() : r.largerEnd();
      b=bb; 
      e=ee;
   }
   //return overlap(r);
   return olp;
}

// both are forward direction
int Range::plusOverlap(const Range &r) const {
	if (e<b || r.e<r.b) return 0;
	int end;
	if (r.b >= b && r.b <= e) {
		end=e<r.e? e : r.e;
		return end-r.b+1;
	}
	else if (b >= r.b && b <= r.e) {
		end = e<r.e? e : r.e;
		return end-b+1;
	}
	else return 0;
}

// both are backward direction
int Range::minusOverlap(const Range &r) const {
	if (b<e || r.b<r.e) return 0;
	Range tmp(e,b);
	return tmp.plusOverlap(Range(r.e,r.b));
}

int Range::plusminusOverlap(const Range &r) const {
	return this->plusOverlap(Range(r.e,r.b));
}
int Range::minusplusOverlap(const Range &r) const {
	Range tmp(e,b);
	return tmp.plusOverlap(r);
}

int Range::overlap(const Range &r, bool sameDirection) const {
	if (*this == r) return length();
	int olp=0;
	if (sameDirection) {
      if (direction() != r.direction()) return 0;
		olp=plusOverlap(r);
		if (olp == 0) {
			olp = minusOverlap(r);
		}
		return olp;
	}
   // only compare opposite directions
	olp = plusminusOverlap(r);
	if (olp == 0) olp=minusplusOverlap(r);
	return olp;
}
int Range::overlay(const Range &r) const {
   int b1=min(b,e);
   int e1=max(b,e);
   int b2=min(r.b,r.e);
   int e2=max(r.b,r.e);
   return min(e1,e2) - max(b1,b2) + 1;
}

bool Range::fuse(const Range &r) {
   if (overlay(r) <= 0) return false;
   int bb=min(smallerEnd(), r.smallerEnd());
   int ee=max(largerEnd(), r.largerEnd());
   b=bb; e=ee;
   return true;
}
void Range::forceFuse(const Range &r) {
   int bb=min(smallerEnd(), r.smallerEnd());
   int ee=max(largerEnd(), r.largerEnd());
   b=bb; e=ee;
}
   
float Range::overlapFraction(const Range &r, bool sd) const {
   // ---------
   // |||||||||
   // ---------
   if (*this == r) return 1;
   //float fr;
   if (sd) {
      if (!sameDirection(r)) return 0;
      if (direction() == '+' && r.direction() == '+') {
         return (min(e,r.e)-max(b,r.b)+1) / (float)(max(e,r.e)-min(b,r.b)+1);
      }
      //fr =  (min(b,r.b)-max(e,r.e)+1) / (float)(max(b,r.b)-min(e,r.e)+1);
      //cerr << *this << " overlaps on - " << r << " by " << fr << endl;
      return (min(b,r.b)-max(e,r.e)+1) / (float)(max(b,r.b)-min(e,r.e)+1);
   }
   // now, only compare Ranges in different directions
   if (sameDirection(r)) return 0;
   int b1, e1, b2, e2;
   b1=min(b,e);
   e1=max(b,e);
   b2=min(r.b, r.e);
   e2=max(r.b, r.e);
   return (min(e1,e2)-max(b1,b2)+1) / (float)(max(e1,e2)-min(b1,b2)+1);
}
float Range::overlayFraction(const Range &r) const {
   int b1, e1, b2, e2;
   b1=min(b,e);
   e1=max(b,e);
   b2=min(r.b, r.e);
   e2=max(r.b, r.e);
   return (min(e1,e2)-max(b1,b2)+1) / (float)(max(e1,e2)-min(b1,b2)+1);
}

// not this is the difference of the two ranges
int Range::distance(const Range& r) const {
   int b1=min(b,e);
   int e1=max(b,e);
   int b2=min(r.b,r.e);
   int e2=max(r.b,r.e);
   return max(b1, b2)-min(e1, e2);
}

// if range contains only one point, then the direction
// is considered unknown
bool Range::sameDirection(const Range &r) const {
   if ( (b<e && r.b < r.e) || (b>e && r.b>r.e)) return true;
   return false;
}
/* 1 HT, 2 HH, 3 TT */
int Range::topology(const Range &r) const {
   if (operator==(r)) return -1;
   if (sameDirection(r)) return 1;
   if (direction() == '+') {
      if (operator<(r)) return 3; // TT
      else if (r < *this) return 2; // HH
      else return 0; // superimpose in opposite direction
   }
   else { // <---  --->
      if (less(r)) return 2;
      else if (r < *this) return 3;
      else return 0;
   }
}

bool Range::plusContain(const Range &r) const {
	if (e<b || r.e<r.b) return false;
	return b<=r.b && e>=r.e;
}
bool Range::minusContain(const Range &r) const {
	if (b<e || r.b<r.e) return false;
	return e<=r.e && b>=r.b;
}
bool Range::plusminusContain(const Range &r) const {
	Range tmp(r.e, r.b);
	return plusContain(tmp);
}
bool Range::minusplusContain(const Range &r) const {
	Range tmp(r.e,r.b);
	return minusContain(tmp);
}

bool Range::contain(const Range &r, bool sameDirection) const {
	bool tmp=(plusContain(r) || minusContain(r));
	if (sameDirection) return tmp;
	else 
		return tmp || plusminusContain(r) || minusplusContain(r);
}

bool Range::contain(int p) const {
   if (direction() == '+') return p >= b && p <= e;
   else return p >= e && p <= b;
}

/* compare with no retard to direciton */
int Range::compare(const Range &r) const {
	if (b==r.b && e==r.e) return 0;
	int b1 = b<e? b : e;
	int e1 = e>b? e : b;
	int b2 = r.b<r.e? r.b : r.e;
	int e2 = r.e>r.b? r.e : r.b;

	if (b1<b2) return -1;
	else if (b1==b2) {
		if (e1<e2) return -1;
		else return 1;
	}
	else return 1;
}

bool Range::lessByDirection(const Range *r) const
{
   if (*this == *r) return false;
   if ((direction() == '+' && r->direction() == '+')
      || (direction() == '+' && r->direction() == '?')
      || (direction() == '?' && r->direction() == '+')
      ) {
      return compare(r) == -1;
   }
   else if (
      (direction() == '-' && r->direction() == '-')
      ||
      (direction() == '?' && r->direction() == '-')
      ||
      (direction() == '-' && r->direction() == '?')) {
      return compare(r) == 1;
   }
   else if (direction() == '+' && r->direction() == '-')
      return true;
   else if (direction() == '-' && r->direction() == '+')
      return false;
   else if (direction() == '?' && r->direction() == '?') {
      cerr << "comparing two points, assuming + direction\n";
      return b<r->b;
   }
   else {
      cerr << "entered impossible state inside lessByDirection()\n";
      exit(1);
   }
}

bool Range::lessByDirection(const Range &r) const
{
   if (*this == r) return false;
   if ((direction() == '+' && r.direction() == '+') 
         || (direction() == '+' && r.direction() == '?') 
         || (direction() == '?' && r.direction() == '+') 
         //|| (direction() == '?' && r.direction() == '?')) 
       )
   {
      return compare(r) == -1;
   }
   else if (
      (direction() == '-' && r.direction() == '-')
      || (direction() == '?' && r.direction() == '-')
      || (direction() == '-' && r.direction() == '?')) {
      return compare(r) == 1;
   }
   else if (direction() == '+' && r.direction() == '-')
      return true;
   else if (direction() == '-' && r.direction() == '+')
      return false;
   else if (direction() == '?' && r.direction() == '?') {
      cerr << "comparing two points inside lessByDirection(). Assuming + direction\n";
      return b<r.b;
   }
   else {
      cerr << "entered impossible state inside lessByDirection(const Range&)\n";
      exit(1);
   }
}

bool Range::lessByDirection(int p, char direc) const {
   if (direction() == '+') {
      return b<p;
   }
   else if (direction() == '-') {
      return b>p;
   }
   else { // ? also a point, assume +, could be bad
      if (direc == '+') return b<p;
      else if (direc == '?') {
         cerr << "point direction cannot be ?\n";
         exit(1);
      }
      else return b>p;
   }
}

/* no more duplication of code
bool Range::greaterByDirection(const Range &r) const {
   if (direction() == r.direction()) {
      if (direction() == '+') return compare(r) == 1;
      else if (direction() == '-') return compare(r) == -1;
      else {
         cerr << "impossible state\n";
         exit(1);
      }
   }
   else if (direction() == '+' && r.direction() == '-')
      return false;
   else if (direction() == '-' && r.direction() == '+')
      return true;
   else {
      cerr << "inside lessByDirection() impossible state\n";
      exit(1);
   }
};

bool Range::greaterByDirection(const Range *r) const {
   if (direction() == r->direction()) {
      if (direction() == '+') return compare(*r) == 1;
      else if (direction() == '-') return compare(*r) == -1;
      else {
         cerr << "impossible state\n";
         exit(1);
      }
   }
   else if (direction() == '+' && r->direction() == '-')
      return false;
   else if (direction() == '-' && r->direction() == '+')
      return true;
   else {
      cerr << "inside lessByDirection() impossible state\n";
      exit(1);
   }
};
*/

/**
 * Situation where this object is before r:
 *
 * b---->e                e<-------b               e<-----b
 *  r.b------>r.e             r.e<------r.b   r.e<------------r.b
 *
 *   b---->e             e<----------b
 * r.b------->r.e          r.e<------r.b
 *
 * This method will mix ranges in different directions
 *
 */
int Range::compareDirectional(const Range &r) const {
   if (b < r.b) return -1;
   if (b > r.b) return 1;
   // now b == r.b, we compare the end
   if (e < r.e) return -1;
   if (e > r.e) return 1;
   return 0;
}

int Range::compareByDirection(const Range &r) const {
   if (*this == r) return 0;
   if (direction() == r.direction()) {
      if (direction() == '+') {
         if (b < r.b) return -1;
         if (b > r.b) return 1;
         if (e < r.e) return -1;
         if (e > r.e) return 1;
         cerr << __FILE__ << ":" << __LINE__ << ": no possible compareByDirection()\n";
         return 0; // should not reach here
         //exit(1);
      }
      else { // - direction
         if (b>r.b) return -1;
         if (b<r.b) return 1;
         // b == r.b
         if (e>r.e) return -1;
         if (e<r.e) return 1;
         cerr << "Comparing " << *this << " with "
            << r << " not possible Range::compareByDirection()\n";
         //exit(1);
         return 0;
      }
   }
   else if (direction() == '+') {
      if (r.direction() == '-') return -1;
      if (r.length() == 1) {
         // this:339440->339427   r:339440->339440
         //cerr << " r is a point inside Range::compareByDirection()\n";
         if (r.b > b) return -1;
         return 1;
      }
      else {
         throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR  r has no direction");
      }
   }
   else if (direction() == '-') {
      if (r.direction() == '+') return 1;
      if (r.length() == 1) {
         // this:37482->37479   r:37482->37482
         //cerr << "comparing " << *this << " with "
         //   << r << " r is a point\n";
         if (r.b>=b) return 1;
         return -1;
      }
      else {
         //cerr << "r has no direction\n";
         throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR  r has no direction");
         //exit(1);
      }
   }
   else if (length() == 1) {
      //cerr << "comparing " << *this << " with " << r << endl;
      //cerr << "this object is a point inside compareByDirection()\n";
      if (r.direction() == '+') {
         if (b <= r.b) return -1;
         return 1;
      }
      else if (r.direction() == '-') { // 3-3  3-2
         if (b >= r.b) return -1;
         return 1;
      }
      else {
         cerr << *this << " compareByDirection() " << r << endl;
         cerr << "entered an impossible state inside compareByDirection()\n";
         return 99;
      }
   }
   return 9;
}

bool Range::before(const Range &r) const {
   if (isPlus() && r.isPlus()) {
      return e < r.b;
   }
   else if (isMinus() && r.isMinus()) {
      return e > r.b;
   }
   else if (isPlus() && r.isMinus())
      return true;
   return false;
}
bool Range::after(const Range &r) const {
   if (isPlus() && r.isPlus()) {
      return b > r.e;
   }
   else if (isMinus() && r.isMinus()) {
      return b < r.e;
   }
   else if (isMinus() && r.isPlus())
      return true;
   return false;
}
   

bool Range::similar(const Range &r, bool sameDirection) const {
	if (sameDirection && direction() != r.direction() )
		return false;

	int error=(int)(marginFraction*(length()<r.length()? length() : r.length()));
	error = error < margin? error : margin;
	if (sameDirection) {
		if (abs(b-r.b)<error && abs(e-r.e)<error) 
			return true;
		else return false;
	}
	int b1=b<e? b : e;
	int e1=e>b? e : b;
	int b2=r.b<r.e? r.b : r.e;
	int e2=r.e>r.b? r.e : r.b;
	if (abs(b1-b2) < error && abs(e1-e2)<error) return true;
	else return false;
}

bool Range::similar(const Range &r, const int marg,
      const float frac, bool sameDirection) const {
	if (sameDirection && direction() != r.direction() )
		return false;

	int error=(int)(frac*(length()<r.length()? length() : r.length()));
	error = error < marg? error : marg;
	if (sameDirection) {
		if (abs(b-r.b)<error && abs(e-r.e)<error) 
			return true;
		else return false;
	}
	int b1=b<e? b : e;
	int e1=e>b? e : b;
	int b2=r.b<r.e? r.b : r.e;
	int e2=r.e>r.b? r.e : r.b;
	if (abs(b1-b2) < error && abs(e1-e2)<error) return true;
	else return false;
}


int Range::length() const { 
   if (e>b) {
      return e-b+1; 
   }
   else {
      return b-e+1;
   }
}

Range Range::copyShrinkBothEnds(int del) const {
   if (length()>2*del) {
      if (direction() == '+') {
         return Range(b+del, e-del);
      }
      else {
         return Range(b-del, e+del);
      }
   }
   else {
      return Range();
   }
}

// range could have negative coordinates!
void Range::growEnd(int icr) {
   if (direction() == '+') {
      e += icr;
   }
   else {
      e -= icr;
   }
}

string Range::asDelimitedString(const char delimiter[]) const {
	ostringstream ous;
	ous << b << delimiter << e;
	return ous.str();
}

string Range::asDelimitedString(const char delimiter) const {
	ostringstream ous;
	ous << b << delimiter << e;
	return ous.str();
}

////////////// RangeChain Class ////////////////////

RangeChain::RangeChain(const RangeChain &ch) 
   : Range(ch), chain(), sumolp(ch.sumolp) 
{
   list<Range*>::const_iterator lit;
   for (lit=ch.chain.begin(); lit != ch.chain.end(); ++lit) {
      chain.push_back(new Range(**lit));
   }
}

RangeChain::~RangeChain() {
	list<Range*>::iterator it;
	for (it=chain.begin(); it != chain.end(); it++) {
		delete (*it);
	}
}

RangeChain::RangeChain(const string &raw) {
   vector<string> r=split(raw, ',');
   //for (int i=0; i<r.size(); i++) {
   for (auto& x : r) { // x should be b-e 123-567 format
      string::size_type s = x.find('-');
      if (s == string::npos) {
         cerr << __FILE__ << ":" << __LINE__ << ":ERROR input string must be in the format b1-e1,b2-e2,...\n";
         throw runtime_error("input string for RangeChain " 
               + raw + " not in format b1-e1,b2-e2,...");
      }
      chain.push_back(new Range(stoi(x.substr(0,s)), stoi(x.substr(s+1))));
   }
}

RangeChain& RangeChain::operator=(const RangeChain &rc) {
   if (this != &rc) {
      list<Range*>::iterator it;
      for (it=chain.begin(); it != chain.end(); it++) {
         delete (*it);
      }
      chain.clear();
      list<Range*>::const_iterator i=rc.chain.begin();
      while (i != rc.chain.end()) {
         chain.push_back(new Range(**i));
         ++i;
      }
      sumolp=rc.sumolp;
   }
   return *this;
}

bool lessRangeByDirectional(const Range *r1, const Range *r2) {
   return r1->lessByDirection(r2);
}

/** this one is doing the merging
 * If overlap with any exon then it will merge the underlying exon; otherwise
 * it will append the Range r after the component Range under consideration.
 * This effect is to sort the component Ranges.
 *
 * There is no sorting action.
 */
void RangeChain::add(const Range &r, bool sd) {
	list<Range*>::iterator it, b, del;
   //char dir=direction();
   chain.sort(lessRangeByDirectional);

	it=chain.begin();
	int mergenum=0;
   //bool inserted=false;
	int olp;
	while (it != chain.end()) {
		if ( (olp=(*it)->overlap(r,sd)) > 0 ) {
			sumolp += olp;
			++mergenum;
			if (mergenum==1) {
				(*it)->merge(r, sd);
				b=it;
				++it;
			}
			else if (mergenum > 1) {
				// merging then delete
				del=it;
				++it;
				(*b)->merge(**del);
				delete *del; 
				chain.erase(del);
			}
			else {
				cerr << "this state is not reachable!\n";
				exit(1);
			}
      }
      else if ((*it)->greaterByDirection(r) ) {
         //cerr << "iterator after r can stop now\n";
         break;
      } 
      else {
         ++it;
      }
	}
	if (mergenum == 0) {
      list<Range*>::iterator it;
      it=find_if(chain.begin(), chain.end(), DRgreater(r));
		chain.insert(it, new Range(r));
	}
}

// length after merge
pair<int,int> RangeChain::length() const {
	list<Range*>::const_iterator it=chain.begin();
	int pluslen=0;
	int minuslen=0;
	while (it != chain.end()) {
		if ((*it)->direction() == '+')
			pluslen += (*it)->length();
		else if ((*it)->direction() == '-')
			minuslen += (*it)->length();
		it++;
	}
	//return make_pair<int,int>(pluslen,minuslen);
	return make_pair(pluslen,minuslen);
}

ostream& operator<<(ostream &ous, const RangeChain &rc) {
   if (rc.numberOfRanges()>0) {
      list<Range*>::const_iterator it;
      for (it=rc.chain.begin(); it != rc.chain.end(); it++) {
         ous << **it << " | ";
      }
   }
   else {
      ous << "NULL";
   }
	return ous;
}


/** using a naive algorithm,
 * cold be improved by indexing in the future if performance
 * is a problem.
 */
pair<int,int> RangeChain::exonOverlap(const RangeChain &rc) const {
   // simple minded M x N algorithm
   list<Range*>::const_iterator i, j;
   i=chain.begin();
   j=rc.chain.begin();
   int olp=0, sumolp=0, numolp=0;
   while (i != chain.end()) {
      j=rc.chain.begin();
      while (j != rc.chain.end() && (olp=(*i)->overlap(**j, true)) == 0) {
         j++;
      }
      if (olp > 0) {
         sumolp += olp;
         ++numolp;
      }
      ++i;
   }
   return make_pair(sumolp, numolp);
}

/** will only carry out this operation if the direction of this object is pure,
 * which means that all the ranges in the chain are in the + or - direction,
 * possibily with a few in the uncertain direction ' '
 */
bool RangeChain::exonMerge(const RangeChain &rc) {
   if (direction() == ' ') {
      cerr << "RangeChain has mixed component Ranges, failed to merger with another Range! Crashed exonMerge() operation\n";
      cerr << *this << endl;
      exit(1);
      //return false;
   }
   if (direction() != rc.direction()) {
      cerr << "exonMerge of RangeChain is only merging chains in the same direction for now\n";
      exit(1);
   }
   list<Range*>::const_iterator j;
   for (j=rc.chain.begin(); j != rc.chain.end(); j++) {
      add(**j, true);
   }
   return true;
}

char RangeChain::direction() const {
   list<Range*>::const_iterator it = chain.begin();
   int plus=0, minus=0, uncertain=0;
   while (it != chain.end()) {
      if ((*it)->direction() == '+') plus++;
      else if ((*it)->direction() == '-') minus++;
      else uncertain++;
      ++it;
   }
   /*
   if (uncertain > 0) {
      cerr << "There are " << uncertain << " component Ranges have uncertain directions while calling direction() of the RangeChain class!\n";
   }
   */
   if (plus > 0 && minus == 0) return '+';
   else if (minus > 0 && plus == 0) return '-';
   else return '?';
}
/** must define these operators to do the sorting */
bool RangeLess(const Range &r1, const Range &r2) {
   return r1.directionalLess(r2);
}
bool RangeLessPtr(const Range *r1, const Range *r2) {
   return r1->directionalLess(*r2);
}
void RangeChain::order() {
   //sort(chain.begin(), chain.end(), mem_fun_ref(&Range::directionalLess));
   //chain.sort(mem_fun(&Range::directionalLess));
   chain.sort(RangeLessPtr);
   //sort(chain.begin(), chain.end(), mem_fun_ref(&Range::directionalLess));
}

int RangeChain::outerLength() const {
	list<Range*>::const_iterator it=chain.begin();
	int min=numeric_limits<int>::max(); //99999999999999999999999999999999;
	int max=0;
	while (it != chain.end()) {
		if ((*it)->smallerEnd() < min) min=(*it)->smallerEnd();
		if ((*it)->largerEnd() > max) max=(*it)->largerEnd();
		it++;
	}
   return max-min+1;
}
Range RangeChain::outerRange() const {
   return Range(chain.front()->begin(), chain.back()->end());
}

int RangeChain::sumLength() const {
   list<Range*>::const_iterator i;
   int len=0;
   for (i=chain.begin(); i != chain.end(); i++) {
      len += (*i)->length();
   }
   return len;
}

pair<int,int> RangeChain::intronLength() const {
   list<Range*>::const_iterator i,j;
   int pluslen=0, minuslen=0;
   i=chain.begin();
   if (direction() == '+') {
      while(1) {
         j = i; ++j;
         if (j == chain.end()) break;
         pluslen += (*j)->begin() - (*i)->end() - 1;
         ++i;
      }
   }
   else if (direction() == '-') {
      while (1) {
         j = i; ++j;
         if (j==chain.end()) break;
         minuslen += (*i)->end() - (*j)->begin() - 1;
         ++i;
      }
   }
   else if (chain.size()>2) {
      while ((*i)->direction() == '+') {
         j=i; ++j;
         if (j == chain.end() || (*j)->direction() == '-') 
            break;
         pluslen += (*j)->begin() - (*i)->end() - 1;
         ++i;
      }
      if (j != chain.end()) {
         i=j; ++i;
         while (i != chain.end()) {
            j=i; ++j;
            if (j == chain.end()) break;
            minuslen += (*i)->end() - (*j)->begin() - 1;
            ++i;
         }
      }
   }
   return make_pair(pluslen,minuslen);
}

int RangeChain::maxIntronLength() const {
   int maxlen=0;
   if (chain.size() == 1) return 0;
   list<Range*>::const_iterator i,j;
   i=chain.begin();
   j=i;

   while (i != chain.end() && (*i)->direction() == '+') {
      j=i; ++j;
      if (j == chain.end() || (*j)->direction() == '-') 
         break;
      if ((*j)->begin() - (*i)->end() - 1 > maxlen)
         maxlen = (*j)->begin() - (*i)->end() - 1;
      ++i;
   }
   if (i == chain.end() || j == chain.end()) return maxlen;
   if (j != chain.end()) {
      i=j; 
   }
   while (i != chain.end()) {
      j=i; ++j;
      if (j == chain.end()) break;
      if ((*i)->end() - (*j)->begin() - 1 > maxlen) 
         maxlen = (*i)->end() - (*j)->begin() - 1;
      ++i;
   }
   return maxlen;
}

ostream& RangeChain::tableRows(ostream &ous, const string &prefix) const {
   list<Range*>::const_iterator i;
   int cnt=1;
   for (i=chain.begin(); i != chain.end(); i++) {
      ous << prefix << "\t" << cnt++ << "\t" << (*i)->begin() << "\t" << (*i)->end() << endl;;
   }
   return ous;
}
} // orpara namespace
