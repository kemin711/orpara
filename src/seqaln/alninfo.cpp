#include "alninfo.h"

namespace orpara {
AlignInfo& AlignInfo::operator=(const AlignInfo &ai) {
   if (this != &ai) {
      score = ai.score;
      identity = ai.identity;
      coverage = ai.coverage;
      readlen = ai.readlen;
   }
   return *this;
}
}

