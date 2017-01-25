#ifndef SEQALN_HELPER_H
#define SEQALN_HELPER_H

// (c) 1997 Kemin Zhou at The Moclecular Sciences Institute

namespace orpara {
// they are defined in C++ std lib
// They should not be used
int max(int i1, int i2, int i3) { 
      return max(max(i1,i2),i3); 
} 
int max(int i1, int i2) { 
      if (i1>i2) return i1; 
         return i2; 
} 
}
#endif
