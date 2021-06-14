#include "genomicinterval.h"

namespace orpara {

bool GenomicInterval::operator<(const GenomicInterval& other) const {
   if (refid < other.refid) 
      return true;
   if (refid > other.refid) return false;
   return Interval::operator<(other);
}

bool GenomicInterval::operator>(const GenomicInterval& other) const {
   if (refid > other.refid) return true;
   if (refid < other.refid) return false;
   return Interval::operator>(other);
}

bool GenomicInterval::operator==(const GenomicInterval& other) const {
   if (refid != other.refid) return false;
   return Interval::operator==(other);
}

} // orpara name space ends
