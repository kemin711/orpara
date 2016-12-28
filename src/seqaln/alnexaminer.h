#ifndef ALNEXAMINER_H
#define ALNEXAMINER_H

#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace orpara {
/**
 * [begin, end) identity
 * Demark a alignment setment.
 */
class Alnseg {
   public:
      /**
       * [start, end) range. End is 1 over the end
       */
      Alnseg(int start, int end, float iden)
         : b(start), e(end), identity(iden) { }
      friend ostream& operator<<(ostream &ous, const Alnseg &as) {
         ous << as.b << "-" << as.e << ":" << as.identity;
         return ous;
      }
      int getSegLength() const { return e-b; }
      int length() const { return e-b; }

      int b,e;
      float identity;
};

/**
 * Give the alignment middle line.
 * This object will figure out the alignments from chimeras.
 */
class Alnexaminer {
   private:
      string middle;
      vector<int> wsum;
      vector<int> deriv;
      vector<Alnseg> segment;
      /**
       * Identical char for displaying sequence alignment.
       */
      static char idenchar;
      /** the window size is not supposed to change very oftern */
      static int w;
      static const int zeroZone=10;
      /**
       * Compute the moving sum of window size w
       */
      void buildMovingSum();
      /**
       * Helper function to convert integer values into
       * derivatives.
       * Doing integer differential.
       */
      void computeDerivative();
      bool allZeros(int i);
      vector<int> findTransitionPoint();
      static float computeIdentity(const string &bar);
      void saveSegment(int bb, int ee);

   public:
      Alnexaminer() : middle(), wsum(), deriv(), segment() { }
      Alnexaminer(const string &m) 
         : middle(m), wsum(middle.size()-w+1),
           deriv(middle.size()-w), segment()
      { }
      /**
       * Find and store the boundary in a vectore.
       * @return the alignment segments.
       */
      const vector<Alnseg>& operator()(const string &bars);
      /**
       * Find boundary of alignments
       * @return true if found segments false otehrwise.
       */
      bool findBoundary();
      /**
       * @return The same data as returned by the operator.
       */
      const vector<Alnseg>& getSegment() const { return segment; }
      ostream& printResult(ostream& ous) const;
      /**
       * @return true if 3 or fewer setments and 
       *   at least one segment >0.99 identical to reference.
       */
      bool isChimera() const;
      /**
       * Print all the segments in one line. EOF not printed.
       */
      ostream& printSegment(ostream& ous) const;
      void debugShow(ostream &ous) const;
};
}
#endif
