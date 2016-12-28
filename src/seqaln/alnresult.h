#ifndef ALNRESULT_H
#define ALNRESULT_H

namespace orpara {
/**
 * Alignment Result for consumption by other programs
 * This object will be returned from the alignment method.
 */
class Alnresult {
   private:
      /** This is the preliminary result
       */
      list<pair<int, int> > alnidx;

      /**
       * Score is dependent on the scoring method used
       * This is derived directly from the algorithm run.
       */
      int score;

      // all the following are derived from the 
      // prilimary result.
      int idencnt;
      int numgaps1, numgaps2;
      int gaplen1, gaplen2;
};
}
#endif
