#ifndef INSERTSORTLIST_H
#define INSERTSORTLIST_H

#include <list>

namespace orpara {

using namespace std;

/**
 * This method is good for partially sorted list
 * that is holding very large objects.
 * If cmp function not provided then use the default
 * less operator.
 * Object must have the move constructor implemented.
 */
template<class T, class C=std::less<T>> void insertSortList(std::list<T>& data, C cmp=std::less<T>()) {
   typename std::list<T>::iterator it, itt;
   it = data.begin();
   itt = it; ++itt;
   // sort first two element
   while (itt != data.end()) {
      if (cmp(*it, *itt)) { // special case *it < *itt, nothing needs to be done
         ++itt; ++it;
      }
      else { // need to find location of itt before it
         // it will have to move left at leat one poistion
         while (it != data.begin() && cmp(*itt, *it)) {
            --it;
         }
         if (cmp(*itt, *it)) { // it point to begin
            data.insert(it, std::move(*itt));
            itt = data.erase(itt);
         }
         else { // move itt to after it
            ++it;
            data.insert(it, std::move(*itt));
            itt = data.erase(itt);
         }
         it = itt;
         --it;
      }
   }
}

}
#endif
