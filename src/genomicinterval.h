#ifndef GENOMICINTERVAL_H
#define GENOMICINTERVAL_H

#include "interval.h"
// expands the interval with reference id

namespace orpara {
/**
 * A simple representation of a region on one genomic DNA.
 */
class GenomicInterval : public Interval {
   public:
      GenomicInterval(int ri, int b, int e) 
         : Interval(b,e), refid(ri)
      {  }
      bool operator<(const GenomicInterval& other) const;
      bool operator>(const GenomicInterval& other) const;
      bool operator==(const GenomicInterval& other) const;
      int getRefid() const {
        return refid;
      } 

   protected:
      /**
       * Integer representation of the reference sequence
       * the caller should know how to interpret the id int
       * meaningful names, such as chr1, chrM, chrX, or chrY
       */
      int refid;
};
}
#endif
