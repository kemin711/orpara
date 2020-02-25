#ifndef RANGE_H
#define RANGE_H

// (c) 2002 Kemin Zhou at orapra.com

#include <iostream>
#include <list>
#include <utility>

using namespace std;

namespace orpara {

/** 
 * @brief a directional segment in one dimentioan space
 *
 * The Range class for sequence alignment.
 * Abstraction of two numbers in one dimentional space.
 * Represent a line.
 *
 * Right now in integer space, can be generalized to a template library.
 *
 * [b, e]
 * 0 or 1-based index should all be accomodated.
 *
 * This is a one-dimensional range.
 * It contains the interval situation where
 * [b,e] contain also the directional information.
 * with b<e as + and b>e as -
 */
class Range {
	public:
      /** This is the default object
       * we used [0,0]
       * to represent the NULL default object
       * A point is actually a Range object.
       * This is an authentic object.
       * Actually this class should also support 
       * negative ranges like [-100, -9]
       * But I use this only for positive integers
       * for dealing with protein or DNA coordinates.
       */
		Range() : b(0), e(0) {}
		Range(int begin, int end) : b(begin), e(end) {} 
      Range(pair<int,int> p) : b(p.first), e(p.second) { }
		Range(const Range &r) : b(r.b), e(r.e) {}
      /** special case, one point is also a Range
       * its direction should be determined by context
       * You can still have direction with one point,
       * in DNA sequence the upper and the lower strands
       * have different meaning in upper A in lower it is T.
       */
      Range(int p) : b(p), e(p) { }
      /** copy the b and e members */
		Range& operator=(const Range &r);
		virtual ~Range() { }
      /** shift to the right operation */
      Range operator+(int c) const { return Range(b+c, e+c); }
      /** shift to the left */
      Range operator-(int c) const { return Range(b-c, e-c); }
      friend Range operator-(int c, const Range &r) { return Range(c-r.b, c-r.e); }
      friend Range operator+(int c, const Range &r) { return Range(c+r.b, c+r.e); }
      Range& operator+=(int c) { b += c; e += c; return *this; }
      Range& operator-=(int c) { b -= c; e -= c; return *this; }

      ///////// information about this object /////////
      /** Note the begin is less than end if in +
       * direction; more than end if in - direction
       */
      int begin() const { return b; }
      int end() const { return e; }
      /** An empty range is one whose begin and end
       * are both at the origin.
       */
      bool empty() const { return b==0 && e==0; }
      /** A NULL range is the same as an empty range.
       * Just an alias
       * [0,0] is considered Null or empty
       */
      virtual bool isNull() const { return b==0 && e==0; }
      virtual void setNull() { b=0; e=0; }
      /** the integer length of the range
       * is is the [max(b,e)-min(b,e)+1]
       * This is the overall length.
       */
		int length() const;
		int smallerEnd() const { return b<e? b : e; }
      /** alias for smallerEnd() */
      int low() const { return b<e? b : e; }
		int largerEnd() const { return e>b? e : b; }
      /** alias for largerEnd() */
      int high() const { return e>b? e : b; }
		/** 
       * @return '+' if forward direction elseif '-' if backward direction
		 *        '?' if b==e
       *   '+' is the plus strand, '-' is the lower strand
       *   in molecular biology.
       *   The direction may also be determined by context.
		 */
		char direction() const { return e>b? '+' : (b>e? '-' : '?'); }
      bool isPlus() const { return e>b; }
      bool isMinus() const { return e<b; }
      /** Special case of range is a point */
      bool isPoint() const { return b==e; }

      //////// Relationship with other Ranges /////////////////
      //////// Comparative information        /////////////////

      /** 
       * calculate the distance between two ranges
       * regardless of direction.
       * If there is no overlap of the two Ranges, then
       * it returns a positive number; otherwise, a negative
       * number is returned.
       * Note: this is the difference.  For intron length
       * between two exon, you need to reduce by one.
       */
      int distance(const Range &r) const;

		/**
       * To restrict the operation to a particular directions
       * use plusOverlap or minusOverlap that will test for 
       * direction first.  
       *
       * @param sameDirection dictates the behavior of 
       *    this function: if set true, then only compute overlap
       *    of Ranges in the same direction. Return 0
       *    if the direction of the two Ranges are different.
       *    If set false, then only compare Ranges in the opposite
       *    direction, return 0 if they are in the same direction.
       *
       * this can be a problem. Be careful.
       *
       * @return positive number that measure the shared length of the
       * two ranges if overlap, 0 otherwise.
		 */
		int overlap(const Range &r, bool sameDirection=true) const;
      /** this method is similar to overlap but has very different
       * behaviors.  It finds the overlap of the two ranges
       * regardless of direction. If no overlap a zero or negative
       * number is returned.
       * @return the amount of length shared regardless of direction.
       *   the number is positive if they overlay.
       *   ---->
       *       ----->  zero if on end of the object is the same
       *               as one end of another object.
       *   negative numnber that is the negative of the 
       *   distance between the two objects.
       */
		int overlay(const Range &r) const;

      /** Range r1 b1 ---------- e1
       *                  |||||
       *  Range r2     b2 ---------- e2
       *
       * return the overlap part divided by the total length
       * after the merge b1 ---- e2
       * (e1-b2+1)/(e2-b1+1), there could be other configurations
       *
       * A zero or negative number indicates no overlap.
       *
       * @param sd dictates the behavior of this function
       *    if set to true then only Ranges in the same direction
       *    will be computed.  If two Ranges are in opposite
       *    directions then return 0. Vise versa.
       *    If set to false, only compare ranges in opposite direction!
       */
      float overlapFraction(const Range &r, bool sd=true) const;
      /** this function cumpute the overlap regardless of direction
       * as compared to the combined Range of this and r.
       */
      float overlayFraction(const Range &r) const;
      /**
       * This method assume that the range is a interval (b<e).
       * Compute overlap if both Ranges are in the + direction,
       * otherwise, return 0.
       */
		int plusOverlap(const Range &r) const;
		int minusOverlap(const Range &r) const;
		int plusminusOverlap(const Range &r) const;
		int minusplusOverlap(const Range &r) const;
		/** 
       * this range contains another
       * Both are intervals; they must be both +.
       * Identical case is special case of containment.
       *
       * @return true if this Range contains r
		 */
		bool plusContain(const Range &r) const;
		/** 
       * both must be in the - direction; otherwise
       * return false.
		 */
		bool minusContain(const Range &r) const;
      /** This Range in the + direction; r in the - direction.
       * They then have a containment relationship.
       */
		bool plusminusContain(const Range &r) const;
      /** This Range in the - direction; the other
       * in the plus direction
       */
		bool minusplusContain(const Range &r) const;

		/** 
       * @param sameDirection If true, then only Ranges in the same directions
       * are compared. Under this situation Range in the opposite directions
       * are considered not contained.  If sameDirection is false, then the
       * direction of the range does not matter. 
		 *
		 * Identical Ranges are special cases of containment.
       *
       * Checks that this object contain r
       *  this object [  [ r ] ]
		 */
		bool contain(const Range &r, bool sameDirection=true) const;
      /** special version for single number
       * [3, 7] conain 3 or 7
       */
      bool contain(int p) const;

      /** this range inside another range,
       * regardless of direction. this is the opposite of contain.
       * Identical is special case of inside.
       * */
      bool inside(const Range &r) const {
         return min(b,e) >= min(r.b, r.e) && max(b,e) <= max(r.b, r.e);
      }
		/** if two ranges have both their ends very close to the
		 * margin defined by the margin static variable then we
		 * consider them to be similar
		 *
		 * The behavior of this function is controlled by the
		 * margin (default 20) and marginFraction (default 0.1)
       * static member variables.
       * @param sameDirection. Default true. If true, then only compare
       * Ranges in the same direction, else (false) the direction of the
       * Range does not matter.  
       *
       *   ----->
       *   <------ would be similar if sameDirection == false.
		 */
		bool similar(const Range &r, bool sameDirection=true) const;

      /** this function give the user more control */
		bool similar(const Range &r, const int marg=20, const float frac=0.1,
            bool sameDirection=true) const;

		/** compare two ranges regardless of direction
       * So two ranges are regarded as in the + direction,
       * then compare the smaller end, then larger end
		 * return -1 if this object is smaller than r
       *        0 if identical
       *        1 if greater
       * This is for ordering ranges. When comparing A & B
       *     --A-->
       * <--B---
       * B < A.  
       * Using SQL statements, such as
       * order by tbegin, tend, A < B
		 */
		int compare(const Range &r) const;
		int compare(const Range *r) const { return compare(*r); }
      /** same as compare returning -1
       * Comparing two Ranges regardless of direction.
       * This object as A, 
       * <---A----    A.less(r) ==> true
       *   ----r-->
       */
      bool less(const Range &r) const { return compare(r) == -1; }
      /** 
       * compare begin, then end
       * return -1 if this object is before r
       *         1 if this object is after r
       *         0 if this object is the same as r
       *         9 for debug, if r is a point
       *         99 if *this is a point
       * if the ranges are in the same direction then
       * compare and compareDirectional should give the
       * same result.
       *
       * If the two ranges are in different directions
       * then the result of this and compare give different
       * results.
       *
       * The effect is sort by begin then end, the object
       * occur to the left < right
       *
       * Note: This comparison is different from compareByDirection
       *       where object on + direction is always less than 
       *       objects on - direction.
       */
		int compareDirectional(const Range &r) const;
      /** The ones in the + direction are before those on -
       * order by direction first, then directionalCompare
       * 
       * @return -1 this before r, 1 this after r, 0 identical
       *
       *  ---> before -->
       *  on - strand
       *  <--A-- after <--B-- of B before A
       */
      int compareByDirection(const Range &r) const;

		bool directionalLess(const Range &r) const { return compareDirectional(r) == -1; }

      /** both range must not be a single point,
       * if at least one of them is a point then return false
       */
      bool sameDirection(const Range &r) const;
      /** @return three posibilities:
        * -1 superimpose in the same direction
        * 1 Head-to-Tail same as sameDirection. 
        * 2 Head-to-Head
        * 3 Tail-to-Tail
        * 0 super impose in different direction
        */
      int topology(const Range &r) const;

      /** comparative operators */
      /** test identity */
		bool operator==(const Range &r) const { return b==r.b && e==r.e; }
		bool operator!=(const Range &r) const { return b!=r.b || e!=r.e; }
		/** 
       * This object before r, use compare()
       * The comparision is directionless.
       * */
		bool operator<(const Range &r) const {
			return compare(r) == -1; }
      /** greater */
		bool operator>(const Range &r) const {
			return compare(r) == 1; }
      /** 
       * If this object direction() == + and r +
       * compare -1 
       * if both are - compare 1 (greater)
       * else + before -
       *
       * (1) both in + direction A reprenset this object
       * -----A---->
       *         ---r-->  A.lessByDirection(r) ==> true
       * (2) both on - direction
       *            <----A------
       *     <--r-----           A.lessByDirection(r) ==> true
       * (3) This object on + and r on -
       *                               -----A---->
       *      <---r---- 
       *      A.lessByDirection(r) ==> true
       * (4) This object on + and r on -
       *                               <-----A----
       *      ---r---->
       *      A.lessByDirection(r) ==> false
       *  Basically you travel from the beginning in
       *  + direction then you make a U turn and travel
       *  backward.
       *
       *  There are total 9 possibilities:
       *   3 (+,-,?) for each range
       *
       *   When both ranges are points, then they are
       *   considered in the + direction.
       *   Usually the user should know the direction in 
       *   such cases, special care should be taken.
       *   This is still a problem that I have not solved yet.
       *
       *  >>>>>>>>>>>>>>>>>+direction+>>>>turn back|
       *  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|
       *
       *  This is specially usefull when comparing exons
       *  structure on the genome.
       */
      bool lessByDirection(const Range &r) const;
      bool lessByDirection(const Range *r) const;
      /** This is a special case.
       *   ----A--->
       *   |p         A.lessByDirection(p) => false
       *
       *  @param direc is the assumed direction of the point
       *         default to +. Cannot be ?.
       */
      bool lessByDirection(int p, char direc='+') const;
      /** the opposite of lessByDirection */
      bool greaterByDirection(const Range &r) const { return !lessByDirection(r); }
      bool greaterByDirection(const Range *r) const { return !lessByDirection(r); }
      bool greaterByDirection(int p) const { return !lessByDirection(p); }
      /** ==A=>  ==B==> A before B on +.
       * <==B==  <==A== A before B on - no overlap
       * No overlap is the key of this function.
       * <==B==  ==A==> A before B by direction + before -
       */
      bool before(const Range &r) const;
      /** The opposite of before 
       * ==B=>  ==A==> A after B on +.
       * <==A==  <==B== A after B on - no overlap
       * <==A==  ==B==> A after B by direction - after +
       * No overlap is the key of this function.
       * */
      bool after(const Range &r) const;

      ////////////// operations on objects ///////////////
		/**
       * Merge this range with r only if they overlap.
       *
       * If sd is true, then only Ranges in the same direction
       * will be merged; otherwise, ranges in the opposite
       * directions are merged.
       *
       * This function uses the overlap function to test
       * overlap, so it is not good.
       * after this operation, the direction of the object
       * is preserved.
		 *
       * @param sd same direction merge or not.
		 * @return merge happened or not
		 */
		bool merge(const Range &r, bool sd=true);
      /** use overlay to check if they are fusable
       * only fuse ranges that have overlay > 0
       * the object now is in the '+' direction after
       * this operation. It does not preserve the
       * original direction.
       */
      bool fuse(const Range &r);
      /** fuse without checking overlap 
       * alwayse produce a range in the '+' direction
      */
      void forceFuse(const Range &r);

      /** 
       * Combine this Range with r to produce a larger range
       * regardless of directions of the two ranges.
       *
       * Implement the combine operation,
       * it will produce a range that is the sum of the two
       * return the overlap if they do, 
       * 0 means no overlap.
       * If combining Ranges of differnet direction
       * the resulting range is always in the + direction.
       * If the two inputs are in the same direction then
       * direction is preserved.
       *
       *  ==A==  ==B== combining A & B produce an envelop
       *  ============
       *
       * @return overlap length of the two input objects
       *    this and r.
       */
      int combine(const Range &r);

      /** The following functions alter the Range itself
       * */
      /** make the Range take new values */
      void set(int be, int en) { b=be; e=en; }
      void set(const Range &r) { b=r.begin(); e=r.end(); }
      void setEnd(int en) { e=en; }
      /** len can be negative number. This will make the end less
       */
      void endUpdate(int len) { e += len; }
      void setBegin(int be) { b=be; }
      /** will set begin = begin + len; len can be negative
       */
      void beginUpdate(int len) { b += len; }
      /** change the direction of the Range. To get a copy use
       * getReverse() 
       * @return a reference to this object.
       */
      Range& reverse() { set(e, b); return *this; }
      //void reverse() { set(e, b); }
      /** @return a range in the opposite direction 
        * This is the same as copyReverse().
        */
      Range getReverse() const { return Range(e,b); }
      /** @return a Range object that is the reverse of this one */
      Range copyReverse() const { return Range(e, b); }
      /** return a smaller object whose both ends are shrinken
       * by the same amoutn del. If length of this object is 
       * less than 2*del the return default object [0,0].
       */
      Range copyShrinkBothEnds(int del) const;
      /** expand the end of the range */
      void growEnd(int icr);

      /////////// Output Functions ///////////

		/** for human to read, has some helping symbols
       * The derived class Noschain use a virtual print 
       * function to make the output operator function
       * properly with derived classes.
       */
		friend ostream& operator<<(ostream &ous, const Range &r);
		// output in a row
		ostream& output(ostream &ous, const char delimiter[]="\t") const { 
			ous << b << delimiter << e; return ous; }
		ostream& output(ostream &ous, char delim='\t') const { 
			ous << b << delim << e; return ous; }
      /** return a string of b,e */
		string asDelimitedString(const char delimiter[]=",") const;
		string asDelimitedString(const char delimiter=',') const;

      static void setMargin(const int mg) { margin=mg; }


	protected:
      /**
       * The beginning and end of the range. if b<e then from
       * left to right (+) direction. If b>e the from right to 
       * left (-) direction
       */
		int b, e;  // begin and end

		/** the following two variables are used in similiar()
       * Default falue for margin is 16 necleotides on
       * intiation of this class. This is a class wide 
       * variable.
       * The default values are good for managing genomic DNA
       * ranges where they can be large. For EST assembly
       * they should be set to zero.
       *
       * This is used in the similar() method.
       * For EST assembly 16 is a reasonable number.
		 */
		static int margin; // default 16
      /** default value 0.1
       */
		static float marginFraction; // default 0.1
};


/** helper functions for sorting and find
 * operation on list<Range*>
 *
 * This is an adaptor to Range.lessByDirection().
 *
 * Sort by + direction then - direction.
 */
bool lessRangeByDirectional(const Range *r1, const Range *r2);

/** a helper for find the first Range greater than 
 * a given element
 * I have tried bind2rn(mem_fun(&Range::greaterByDirection)) and it
 * is not working well for the compiler!
 * So I have to write something more.
 *
 * This is a simple function object.
 */
class DRgreater {
   const Range *given;
   public:
      DRgreater(const Range &r) : given(&r) { }
      /** adaptor for Range.greaterByDirection() */
      bool operator() (const Range *r2) const {
         return r2->greaterByDirection(given);
      }
};


/**
 * this class merges overlapping Ranges
 * A non-merged chain can be implemented
 * this class is a merged chain of ranges.
 * If ranges overlap then they are merged into one.
 *
 * Derive this class from the base class Range,
 * the overal range of the chain is a Range
 * [B, E]
 *
 * If you know that your chain components will not overlap
 * each other you should use some other simple data structure
 * such as vector of Range.
 * Right now this class is useful in checking Chimeras.
 */
class RangeChain : public Range {
	public:
		RangeChain() : chain(), sumolp(0) {}
      /** copy constructor of the same type */
		RangeChain(const RangeChain &ch);
      /**
       * This constructor takes a string input.
       * It will parse the strin as a list of Ranges.
       * It is not sorted in any way, but could be in
       * one of the two ways:
       *
       * if on +, then from small to large.
       * if on -, then from large to small.
       * The following is an example of +
       *
       * 342525-342791,342797-343540,343546-345213,348262-349290
       *
       * - direction from large to small
       *
       * 258231-257542,257149-256862,256375-255698
       */
      RangeChain(const string &raw);
		RangeChain(const Range &r) : sumolp(0) { 
			chain.push_back(new Range(r)); }
      /** 
       * deallocate memory in each list node 
       */
		~RangeChain();
      RangeChain& operator=(const RangeChain &rc);

		/** This is the major accumulator function to 
       * build the chain from input Ranges.
       *
       * @param sd controls the behavior of this function. The default
       *    is true. This parameter has the same meaning as 
       *    the one in overlap()
       *
       * This operation will merge the added fragment with
       * existing Ranges in the chain. If sd is true, and
       * the direction of the range element of the chain
       * is in opposite direction as the input range, then
       * this method will simply append the range to the
       * end of the chain. So the chain may end up as 
       * heterogenous. The merging will also have a side
       * effect of sorting if sd is set to true.
       *
       * Before adding each new Range, we dorted the chain
       * with lessByDirection() login.
       * So the chain elements will be -----> then <-----
		 */
		void add(const Range &r, bool sd=true);
      /** return the sum of the overlaps while
       * building this chain. The add function
       * keeps track of the overlaps during the operation.
       * If not using the add function, then this number
       * should be zero.
       */
		int getOverlap() const { return sumolp; }
      /** merges the underlying Ranges if they are in the same direction.
       * @return true if merge success; otherwise false
       */
      bool exonMerge(const RangeChain &rc);

		/** obtain the sum length of all Ranges
		 * the firt length is from the plus direction
		 * the second length is from the minus direction
		 * The length is the length after the merger of
		 * individual ranges.
       *
       * @return length(+,-) directions.
       *
       * The object could hold object in opposite directions.
       * if all the composing nodes are in the + direction
       * then the length for the second component will be 
       * zero.
       *
       * This method should be removed. Replaced with
       * exonLength. This method should be reserved for the
       * parent length returning the overal length.
		 */
		pair<int,int> length() const;
      /** return the overall length of the range
       * This is the original lengh of the parent class.
       */
      int outerLength() const;
      /** return the first and the last as a Range
       */
      Range outerRange() const;
      /** length of all Ranges regardless of direction
       * as oppose to length that returns two lengths
       */
      int sumLength() const;
      /** alias for sumLength()
       */
      int exonLength() const { return sumLength(); }
      /** return the sum of the length of the gaps between
       * all ranges for both the plus and the minus
       * direction.
       */
      pair<int,int> intronLength() const;
      /** maximum intron length */
      int maxIntronLength() const;
      /** @return the outer Range begin */
      int begin() const { return chain.front()->begin(); }
      /** @return the outer Range end */
      int end() const { return chain.back()->end(); }
      /** output using ' | ' as delimiter
       */
		friend ostream& operator<<(ostream &ous, const RangeChain &rc);
      /** compare this object with another RangeChain object
       * and compute the ammount of the overlap of underlying
       * Ranges.
       *
       * @return (first,second) the first is the overlap of 
       *     the component Ranges.  The second value is the 
       *     number of ranges that overlap.
       */
      pair<int,int> exonOverlap(const RangeChain &rc) const;
      /** return the number of Ranges in this chain. */
      int numberOfRanges() const { return chain.size(); }
      /** if all of the component Ranges are in the +
       * direction then return +. If all in the - direction
       * then return -.  If mixed then return ' '
       *
       * This also overwrites the parent direction method.
       */
      char direction() const;
      /** order the range first by direction then by directional range.
       * This is hard to implement.  Need a special function object.
       */
      void order();
      /** produce table rows with prefix as prefix.
       * For example the prefix could be a string 
       * "tdb\tgenomic" this function will produce
       * tdb \t genomic \t exnum \t ex1b \t ex1b
       * tdb \t genomic \t exnum \t ex2b \t ex2b
       * ....
       * tdb \t genomic \t exnum \t exnb \t exnb
       */
      ostream& tableRows(ostream &ous, const string &prefix) const;

      typedef list<Range*>::const_iterator const_iterator;
      /** for iterator throught the chain */
      const_iterator itbegin() const { return chain.begin(); }
      const_iterator itend() const { return chain.end(); }

	protected:
      /** Use the list to kepp track of individual Ranges.
       *
       * Could be from large to small or from small to large.
       * Depends on the direction of the underlying ranges.
       */
		list<Range*> chain;
      /** 
       * this is the sum of overlaps during the add() operation.
       */
		int sumolp; 
};
}
#endif
