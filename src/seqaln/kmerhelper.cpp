#include "kmerhelper.h"
#include <bitset> // for debug print out
#include <cmath>
#include <iterator>
#include <set>
//#include <cassert>
//#include <fstream>
//#include <stdexcept>

namespace orpara {
int b2i(char ch) {
   switch (ch) {
      case 'A':
         return 0;
      case 'a':
         return 0;
      case 'C':
         return 1;
      case 'c':
         return 1;
      case 'G':
         return 2;
      case 'g':
         return 2;
      case 'T':
         return 3;
      case 't':
         return 3;
      default:
         //cerr << "base " << ch << " is not ACGT\n";
         //throw range_error("base " + ch + " ignored in kmer");
         return -1;
   }
}

unsigned int computeMask(int mersz) {
   return (2<<(2*mersz-1)) - 1;
}

/**
 * For external use
 */
vector<int> hashArray(const string &s, const int w) {
   unsigned int v=0;
   vector<int> tmp(s.length()-w+1);
   // mask is 00011-w-111 total 64 bits
   unsigned int M=(2<<(2*w-1))-1;
   size_t i;
   for (i=0; i<w; ++i) {
      v <<=2; 
      v |= b2i(s[i]);
   }
   cout << "mask value: " << bitset<32>(M) << endl;
   for (i=0; i<s.length()-w; ++i) {
      //cout << bitset<32>(v) << endl;
      tmp[i]=v;
      v <<=2;
      v |= b2i(s[i+w]);
      v &= M;
   }
   tmp[i]=v;
   cout << "length of string: " << s.length()
      << " length of hash: " << tmp.size() << endl;
   return tmp;
}

// used in kmercount
double computeD2(const vector<double> &f1, const vector<double> &f2) {
   double d=0;
   for (size_t i=0; i<f1.size(); ++i) {
      d += pow(f1[i]-f2[i], 2);
   }
   return sqrt(d);
}

void trimLonely(set<int> &loc) {
   cout << "sorted single location\n";
   copy(loc.begin(), loc.end(), ostream_iterator<int>(cout, ", "));
   cout << endl;
   // remove lone one, that are d away from other points
   auto it=loc.begin();
   set<int>::iterator itt, del, it3;
   static const int d = 10;
   while (it != loc.end()) {
      itt = it; ++itt;
      if (itt == loc.end()) break;
      if ((*it) < (*itt)-d) {
         del=it;
         ++it;
         loc.erase(del);
         continue;
      }
      else break;
   }

   while (it != loc.end()) {
      itt = it; ++itt;
      if (itt == loc.end()) break;
      it3 = itt; ++it3;
      if (it3 == loc.end()) break;
      //cout << "evaluating " << (*it) << " | " << (*itt) << " | " << (*it3) << endl;
      if ((*itt) > (*it)+d && (*itt) < (*it3)-d) {
         cout << "eliminating " << (*itt) << endl;
         loc.erase(itt);
      }
      ++it;
   }
   cout << "\nremoved lonely single location\n";
   copy(loc.begin(), loc.end(), ostream_iterator<int>(cout, ", "));
   cout << endl;
}

void check_pair(set<int> &loc){

	set<int>::iterator it;
     
          

        it=loc.begin();
         while (it != loc.end()) {

     		cout<<*it<<endl;
                ++it;

         }


}


// with this, you don't need to remove
// lonely points
vector<pair<int,int> > combineRanges(set<int> &loc, int d) {
   
   set<int>::iterator it, it2, it3;
   it = loc.begin();
   pair<int, int> res;
   vector<pair<int,int> > largeRanges; // merged regions
   while (it != loc.end()) {

    //  cout<<*it<<endl;
      res.first = *it;
      ++it;
      if (it == loc.end()) break; // end of input
	//cout<<*it<<endl;

      if ((*it) - res.first < d) {
         res.second = *it;
      }
      else {
         ++it; continue;
      }
      ++it;
     
      while (it != loc.end() && (*it) - res.second < d) {
         
         res.second = *it;
         ++it;
      //   cout<<*it<<endl;

      }
      //if (res.second - res.first > 0) {
         //cout << res.first << "-" << res.second << endl;
      largeRanges.push_back(res);
      //}
   }
   return largeRanges;
}


vector<pair<int,int> > combineStemRanges(set<int> &loc) {
   
   set<int>::iterator it, it2, it3;
   it = loc.begin();
   pair<int, int> res;

   int d=1;
   
   vector<pair<int,int> > largeRanges; // merged regions
   while (it != loc.end()) {

    
      res.first = *it;
      ++it;
      if (it == loc.end()) break; // end of input
	
      if ((*it) - res.first < d) {
         res.second = *it;
      }
      else {
         ++it; continue;
      }
      ++it;
     
      while (it != loc.end() && (*it) - res.second < d) {
         
         res.second = *it;
         ++it;
     
      }
     
      largeRanges.push_back(res);
      
   }
   return largeRanges;
}



pair<int,int> longestRegion(const vector<pair<int,int> > &reg) {
   int len=reg.front().second - reg.front().first;
   pair<int, int> maxr=reg.front();
   for (size_t i=1; i<reg.size(); ++i) {
      if (reg[i].second - reg[i].first > len) {
         len = reg[i].second - reg[i].first;
         maxr = reg[i];
      }
   }
   return maxr;
}

pair<int,int> nearestRegion(const vector<pair<int,int> > &reg, int point) {
   int dist=9999999;
   pair<int,int> bestr= {9999, 99999};
   for (size_t i=0; i<reg.size(); ++i) {
      if (point <= reg[i].second && point >= reg[i].first) {
         dist=0;
         return reg[i];
      }
      if (point > reg[i].second) {
         if (point-reg[i].second < dist) {
            dist = point - reg[i].second;
            bestr = reg[i];
         }
      }
      else if (point < reg[i].first) {
         if (reg[i].first - point < dist) {
            dist = reg[i].first - point;
            bestr = reg[i];
         }
      }
   }
   return bestr;
}
}



