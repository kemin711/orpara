#ifndef INTERVAL_H
#define INTERVAL_H

// (c) 2002 Kemin Zhou at orapra.com

#include <iostream>
#include <list>
#include <vector>
#include <cstdlib>
#include <utility>

using namespace std;

namespace orpara {

/** 
 * @brief Reprsenting all foward ranges in one dimentioanl space.
 *
 * the Interval class is a lot simpler than the Range class
 * that considers all 4 possibility between two ranges
 * when comparing them.  In reality, most work is done 
 * with the ranges are int the + direction such as protein
 * alignments. This class will make things simple and fast.
 * b<e
 *
 * Represets [b, e] closed interval where b <= e
 */
class Interval {
   public:
      /**
       * Default constructor.
       * The region is a single point (-1,-1)
       * which is most unlikely number. Could consider
       * min(int) smallerst number.
       * isNull() will test against this value.
       * Maybe should use -inf,-inf for NULL
       */
      Interval() : b(-1), e(-1) { }
      /**
       * Constructor from two numbers.
       * @param bb begin of interfal, smaller end
       * @param ee end of interval, larger end
       */
      Interval(const int bb, const int ee) : b(bb), e(ee) { }
      /**
       * Copy constructor
       */
      Interval(const Interval &iv) : b(iv.b), e(iv.e) { }
      /** use string as input to construct a interval
       * You can use a,b, a-b, a:b formats
       */
      Interval(const string &str);
      /**
       * Constructor from a pair of values.
       */
      Interval(const pair<int,int> &p) : b(p.first), e(p.second) { }
      /**
       * Destructor
       */
      virtual ~Interval() { }
      /**
       * assignment operator
       */
      Interval& operator=(const Interval& other);

      /**
       * Output stream operator
       */
      friend ostream& operator<<(ostream& ous, const Interval& iv) {
         ous << iv.b << "-" << iv.e; return ous; }

      /**
       * If this object is null, then this function
       * return the length of i. So null objects
       * always overlap with other objects. This is
       * useful for looping start conditions.
       *
       * The reverse situation, i is null then  return
       * the length of this object.  If both objects
       * are null then return 0;
       *
       * @return a number representing the length of the 
       *   overlap between the two intervals.
       *   If there is no overlap, then a negative number will be 
       *   returned to represent the distance between the intervals.
       */
      int overlap(const Interval &i) const;
      int overlap(const pair<int,int> &i) const {
         return overlap(i.first, i.second);
      }
      int overlap(const int bb, const int ee) const;
      bool contains(const Interval& other) const {
         if (other.isNull()) return true;
         return b <= other.b && e >= other.e;
      }
      /**
       * @return true if this interval contains the point
       */
      bool contains(int point) const {
         return point >= b && point <= e;
      }
      /**
       * the distance between the two ranges
       *    |Range1|      | Range 2 |
       *           |<-d-->|
       * if overlap then return a negative number as the
       * amount of overlap.
       * if one of the object is null, then return 0.
       */
      int distance(const Interval& i) const;
      /** 
       * Merging two Intervals.
       * The merge only happens when the overlap > cut
       * If this object is null, then the overlap will
       * be the length of i and this object will become i.
       * @param cut the cutoff value for merging.
       * @return the overlapping value regardless of the merging 
       *    happened or not.
       *
       * It is better to use the overlap method to control
       * the merge method.
       */
      virtual int merge(Interval *i, const int cut=5);
      /**
       * Merge this interval with i; this object will
       * be extended to include both the original and the
       * input interval. If this object is Null then
       * this object will become i.
       * @return the overlap value before the two objects merged. 
       *    Negative indicates that there is no overlap.
       */
      int extend(const Interval &i);
      virtual string toDelimitedString(const string &dl="\t") const;
      void clear() {
         b=-1; e=-1;
      }
      /**
       * This object will be come an empty object.
       */
      void setNull() {
         b=-1; e=-1;
      }
      /**
       * This interval should be arranged before iv
       * in sorted order.
       * @return true if this object is before iv.
       */
      bool less(const Interval &iv) const;
      /**
       * Less operator.
       */
      bool operator<(const Interval &iv) const { return less(iv); }
      bool operator>(const Interval &iv) const { 
         if (b>iv.b) return true;
         if (b<iv.b) return false;
         return e > iv.e;
      }
      /**
       * For comparison with single point
       */
      bool operator<(int p) const {
         return e < p;
      }
      friend bool operator<(const Interval& iv, int p) {
         return p < iv.b;
      }
      /**
       * Equal operator.
       */
      bool operator==(const Interval &iv) const { return b==iv.b && e==iv.e; }
      bool operator!=(const Interval &iv) const { return b!=iv.b || e!=iv.e; }
      bool before(const Interval &iv) const { return e < iv.b; }
      bool after(const Interval &iv) const { return b > iv.e; }
      /**
       * @return true if the range is [-1, -1]. 
       *   This is the original default state.
       */
      bool isNull() const { 
         return b==-1 && e==-1; 
      }

      void setBegin(const int bb) { b=bb; }
      void setEnd(const int ee) { e=ee; }
      void set(const Interval &iv) { b=iv.b; e=iv.e; }
      /**
       * change the interval to be [bb, ee]
       */
      void set(const int bb, const int ee) { b=bb; e=ee; }
      /** return the beginning of the interval.
       * begin() will do the same */
      int getBegin() const { return b; }
      /** returns the start of the interval */
      int begin() const { return b; }
      /**
       * @return end of the interval
       */
      int getEnd() const { return e; }
      /**
       * @return end of the interval
       */
      int end() const { return e; }

      /** 
       * more meaningful for composite classes derived from
       * this base class. Such as IntervalChain and IntervalPile
       * @return the number of elements inside this element.
       */
      virtual int size() const { return 0; }
      /** return the length of the Interval
       * For composite derived class, it is the Interval
       * covered by all the members.
       * it is min(b) - max(e) of all members.
       * @return length of the interval. For NULL ranges
       *   the length is defined as zero not 1.
       */
      int length() const { if (isNull()) return 0; return e-b+1; }

   private:
      int b, e;
};

/** 
 * Represents all intervals that pile up on top of each other.
 * It is used to give an average information of a cluster of 
 * overlapping intervals.
 *
 * This class is a container class.
 *
 * The outer Range will contain all component Intervals.
 */
class IntervalPile : public Interval {
   public:
      IntervalPile() : Interval(), members() { }
      IntervalPile(Interval* ip) : Interval(*ip), members() { 
           members.push_back(ip); }
      /** destroies each member, deallocate memory for each
       */
      ~IntervalPile();
      /** merging two piles.
       * This method will combine the members from ip.
       * It will not make new copies because it does not know
       * the type of the under members. It also does
       * a simple transfer of ownership, which means that
       * the container in ip will be empty after this
       * operation.  ip will become 0 after the operation.
       *
       * @param ip must be the Interval& type for the
       * compiler to recognize it as a virtual function!
       *
       * after function call ip becomes 0.
       *
       * No merge operation at element level, just
       * combining the members in two sets.
       */
      int merge(Interval *ip, const int cut=5);
      /** 
       * Add one more interval to this pile.
       * The Interval must be the same type as the
       * members. We are using the base class pointer
       * for polymorphism behavior.
       *
       * @param i input Interval. After function call i becomes 0.
       */
      void add(Interval *i);
      /** the number of Intervals in this object
       */
      int size() const { return members.size(); }
      string toDelimitedString(const string &dl="\t") const;

   protected:
      /** this variable holds all the Interval objects
       * and its derived classes. The derived classes
       * can add extra information to the Interval.
       */
      vector<Interval*> members;
};

/** 
 * Class HalfAlignInterval to represent half of the alignment.
 * This class is used for chimera detection.
 *
 * It extends Interval class.
 *
 * It add identical, alnlen, and cov information
 * for each alignment.
 */

class HalfAlignInterval : public Interval {
   public:
      HalfAlignInterval(const int bb, const int ee, int iden, int alen, float cc)
         : Interval(bb,ee), identical(iden), alnlen(alen), cov(cc) { }
      HalfAlignInterval(const HalfAlignInterval &i) 
         : Interval(i), identical(i.identical), alnlen(i.alnlen), cov(i.cov) { }
      int getIdentical() const { return identical; }
      float getIdentity() const { return static_cast<float>(identical)/alnlen; }
      int getAlnlen() const { return alnlen; }
      float getCov() const { return cov; }

   private:
      /** number of identical matches */
      int identical; 
      /** length of the alignment */
      int alnlen;
      /** How much this alignment covers the covers the 
       * sequence in question. A sequence with length >= alnlen
       * is implied by this class. 
       *
       * This should be the target coverage
       * */
      float cov;
};

/** This is more at application level class.
 * Represents all alignments clusterd together.
 */
class HalfAlignPile : public IntervalPile {
   public:
      HalfAlignPile() : IntervalPile(), identity(-1.0), qcov(-1.0),
         tcov(-1.0) { }
      HalfAlignPile(const HalfAlignPile &p) 
         : IntervalPile(p), identity(p.identity), qcov(p.qcov) { }
      HalfAlignPile(HalfAlignInterval* hp) 
         : IntervalPile(hp), identity(-1.0), qcov(-1.0), tcov(-1.0) { }
      /** return the average identity over the covered regions
       * It is simply computed by sum(identical)/sum(alnlen)
       */
      float getIdentity() const;
      /** return the average target coverage for this pile's 
       *
       * Simple average of all, with respect to target sequence.
       */
      float getCov() const { return getTcov(); }
      /** computed by the sum of all intervals divided
       * by the outer most interval length of this pile
       * divided again by the number of intervals in 
       * this pile
       * This is called the local_qcov in the context 
       * of the chain.
       */
      float getQcov() const;
      float getTcov() const;
      /** print the information in tabular format.
       * The columns:
       * aln_begin, aln_end, pile_size, avg_identity, 
       *    local_qcov, avg_tcov
       * 
       * @param dl is the field separator
       */
      string toDelimitedString(const string &dl="\t") const;

   private:
      /** helper function used by this object 
       * to compute identity, qcov, and tcov
       * It only changes mutable variables and looks
       * like a const function.
       */
      void computeIntermediateResult() const;
      /** for repeated calls of getting the average identity 
       * This is a tmp variable for getIdentity()
       * **/
      mutable float identity;
      /** tmp variable for getCov()
       * Average query coverage on this pile 
       **/
      mutable float qcov;
      /** average tcov **/
      mutable float tcov;
};

/*

bool lessInterval(const Interval *iv1, const Interval *iv2) {
   return iv1->less(*iv2);
}
*/

/** 
 * IntervalChain: A composition of Interval.
 * This class used a list of Interval as member to store
 * succeesive alignment information on one half of the
 * matches.  It is intended to be used for chimera detection
 * of a large alignment input table.
 *
 * Derived class of Interval if implementing the merge method
 * can use this class to manage more information.  For example,
 * if you want to count the number of merged alignment,
 * the average coverage and identity you can make a derived
 * class of the Interval.
 */

class IntervalChain : public Interval {
   public:
      typedef list<Interval*>::iterator iterator;
      typedef list<Interval*>::const_iterator const_iterator;
      iterator begin() { return chain.begin(); }
      iterator end() { return chain.end(); }
      IntervalChain() : Interval(-1, -1), chain() { }

      IntervalChain(const int bb, const int ee) 
         : Interval(bb, ee), chain() { chain.push_back(new Interval(bb,ee));  }
      IntervalChain(const Interval &iv) : Interval(iv), chain()
         { chain.push_back(new Interval(iv)); }
      IntervalChain(Interval* iv) : Interval(*iv), chain() { chain.push_back(iv); }
      ~IntervalChain();

      /** the add operation keeps the chain sorted
       * It also update the overal range.
       * This method should work with derived classes.
       * For polymorphic behavior.
       * This function will make a copy of *iv
       *
       * After the operation, iv becomes 0
       *
       * The object pointed to by iv is the same type
       * as those in the chain.
       * Use the merge operation to merge overlapping
       * intervals.
       */
      void add(Interval *iv);
      /** this method will work with the bare minimus for speed. 
       * This methods is used when you only care about the 
       * range with no auxiliary information such as identity, etc.
       */
      void add(const int bb, const int ee);
      int size() const { return chain.size(); }
      /* output function */
      friend ostream& operator<<(ostream &ous, const IntervalChain &iv);

      /** sume of the qcov of each pile 
       * sum(pile length), seq len stored outside of this object. 
       * */
      int getInnerLength() const;
      /** outer Range of the chain: min to max */
      int getOuterLength() const;

   private:
      /** keep a sorted set of non-mergeble Intervals
       * Sorted list, the merge operation will do the sorting
       * implicitly.
       * The list is usually very small 1 to 3 at most
       * 99% size=1
       */
      //set<Interval*, lessInterval> chain;
      list<Interval*> chain;
};
}
#endif
