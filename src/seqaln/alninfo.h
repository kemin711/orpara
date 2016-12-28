#ifndef ALIGNINFO_H
#define ALIGNINFO_H

#include <iostream>

using namespace std;

namespace orpara {
/**
 * Simple number for the alighment.
 */
class AlignInfo {
   public:
      AlignInfo(int sc, float iden, float cov, int rawlen)
         : score(sc), identity(iden), coverage(cov), readlen(rawlen) { }
      AlignInfo(const AlignInfo &ai) 
         : score(ai.score), identity(ai.identity), coverage(ai.coverage),
            readlen(ai.readlen) { }
      AlignInfo& operator=(const AlignInfo &ai);
      bool operator<(const AlignInfo &ai) const {
         if (score < ai.score) return true;
         if (score > ai.score) return false;
         if (identity < ai.identity) return true;
         if (identity > ai.identity) return false;
         return readlen < ai.readlen;
      }

      bool operator==(const AlignInfo &ai) const {
         return score == ai.score && identity == ai.identity
            && coverage == ai.coverage && readlen == ai.readlen;
      }
      friend ostream& operator<<(ostream &ous, const AlignInfo &ai) {
         ous << ai.score << '\t' << ai.identity << '\t' << ai.coverage
            << '\t' << ai.readlen;
         return ous;
      }

   private:
      int score;
      float identity;
      float coverage; // coverage of the ref
      int readlen;
};
}
#endif
