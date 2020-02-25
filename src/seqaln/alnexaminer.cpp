#include "alnexaminer.h"
#include <iterator>
#include <stdexcept>

namespace orpara {
char Alnexaminer::idenchar='|';
int Alnexaminer::w=30;

void Alnexaminer::buildMovingSum() {
   int movingSum=0;
   int i;
   for (i=0; i < w; ++i) {
      if (middle[i] == idenchar) ++movingSum;
   }
   i=0;
   while (i < (int)(middle.length()-w)) {
      wsum[i]=movingSum;
      if (middle[i+w] == idenchar) ++movingSum;
      else --movingSum;
      if (middle[i+1] == idenchar) --movingSum;
      else ++movingSum;
      ++i;
   }
}

bool Alnexaminer::allZeros(int i) {
   int leftIdx=(i-zeroZone) < 0 ? 0 : i-zeroZone;
   int rightIdx=(i+zeroZone) > (int)deriv.size() ? deriv.size() : i+zeroZone;
   int k;
   for (k=leftIdx; k<i; ++k) {
      if (deriv[k] != 0) {
         return false;
      }
   }
   // the left are all zeros in zero zone
   for (k=i+1; k < rightIdx; ++k) {
      if (deriv[k] != 0) {
         return false;
      }
   }
   return true;
}

vector<int> Alnexaminer::findTransitionPoint() {
   vector<int> transition;
   int i=0;
   while (i < int(deriv.size()-w)) {
      if (deriv[i] != 0) {
         if (allZeros(i)) { // not a boundary point
            i += zeroZone; continue;
         }
         int sumL=0;
         int sumR=0;
         int Ileft = (i-w < 0) ? 0 : i-w;
         int Iright = (i+w) > (int)deriv.size() ? deriv.size() : i+w;
         for (int j=Ileft; j<i; ++j) sumL += deriv[j]*deriv[j];
         for (int j=i+1; j<Iright; ++j) sumR += deriv[j]*deriv[j];
         if (abs(sumL - sumR) > 1) {
            transition.push_back(i);
            i += w;
         }
         else ++i;
      }
      else ++i;
   }
   return transition;
}

void Alnexaminer::computeDerivative() {
   for (size_t i=0; i<wsum.size()-2; ++i) {
      deriv[i]=wsum[i+1]-wsum[i];
   }
}

float Alnexaminer::computeIdentity(const string &bar) {
   int S=0;
   for (size_t i=0; i<bar.length(); ++i) {
      if (bar[i] == idenchar) ++S;
   }
   return float(S)/bar.length();
}

int pairDiff(const pair<int,int> &p) {
   return p.second-p.first;
}

void Alnexaminer::saveSegment(int bb, int ee) {
   if (bb > (int)middle.length()) {
      cerr << "bb outof range\n";
   }
   segment.push_back(Alnseg(bb, ee, computeIdentity(middle.substr(bb, ee-bb))));
}

bool Alnexaminer::findBoundary() {
   buildMovingSum();
   computeDerivative();
   vector<int> transit=findTransitionPoint();
   //copy(transit.begin(), transit.end(), ostream_iterator<int>(cout, ", "));
   //cout << endl;
   // collect raw results
   vector<pair<int,int> > range;
   int i=0;
   while (i < (int)transit.size()) {
      int start=transit[i];
      int e=start;
      ++i;
      while (i < (int)transit.size() && transit[i] - e < 3*w) {
         e = transit[i]; ++i;
      }
      range.push_back(make_pair(start,e));
   }
   if (range.empty()) {  // some debug output
      return false;
   }
   // move result to the nearest segment boundary
   for (i=0; i < (int)range.size(); ++i) {
      range[i].first += w;
      while (range[i].first < (int)middle.length() && middle[range[i].first] != ' ') 
         ++(range[i].first);
      range[i].second += w;
      if (range[i].second > (int)middle.length()) {
         range[i].second=middle.length();
      }
      if (int(middle.size())-int(range[i].second) < 3*w) {
         range[i].second = middle.size();
      }
      else {
         range[i].second += 2*w;
         if (range[i].second > (int)middle.length()) 
            range[i].second = middle.length();
         while (middle[range[i].second] != ' ') --range[i].second;
         ++(range[i].second);
      }
   }
   // remove very closely packed regions, cutoff w
   if (range.size() > 1) {
      vector<pair<int,int> > combined(1, range[0]);
      i=1;
      while (i < (int)range.size()) {
         if (range[i].first - combined.back().second < w) {
            combined.back().second = range[i].second;
         }
         else {
            combined.push_back(range[i]);
         }
         ++i;
      }
      if (range.size() != combined.size())
         range=combined;
   }
   // transform into final results
   if (range[0].first > 0 && range[0].first < 2*w) { 
      range[0].first=0; 
   }
   int lastEnd=0;
   for (i=0; i < (int)range.size(); ++i) {
      if (range[i].first - lastEnd > 0) {
         saveSegment(lastEnd, range[i].first);
      }
      saveSegment(range[i].first, range[i].second);
      lastEnd=range[i].second;
   }
   if (lastEnd < (int)middle.length()) {
      saveSegment(lastEnd, middle.length());
   }
   return true;
}

const vector<Alnseg>& Alnexaminer::operator()(const string &bars) {
   middle=bars;
   wsum.resize(middle.size()-w+1);
   deriv.resize(middle.size()-w);
   segment.clear();
   findBoundary();
   return segment;
}

ostream& Alnexaminer::printResult(ostream& ous) const {
   for (size_t i=0; i<segment.size(); ++i) {
      ous << segment[i] << endl
         << middle.substr(segment[i].b, segment[i].getSegLength())
         << "*END\n";
   }
   return ous;
}

ostream& Alnexaminer::printSegment(ostream& ous) const {
   ous << segment.front();
   for (size_t i=1; i<segment.size(); ++i) {
      ous << '\t' << segment[i];
   }
   return ous;
}

bool Alnexaminer::isChimera() const {
   if (segment.empty() || segment.size()>3) return false;
   for (size_t i=0; i<segment.size(); ++i) {
      if (segment[i].length() > 200 && segment[i].identity>0.99) 
         return true;
   }
   return false;
}

void Alnexaminer::debugShow(ostream &ous) const {
   ous << "middle line from alignment:\n"
      << middle << endl
      << "moving window sum:\n";
   copy(wsum.begin(), wsum.end(), ostream_iterator<int>(ous, " "));
   ous << endl << "derivative:\n";
   copy(deriv.begin(), deriv.end(), ostream_iterator<int>(ous, " "));
   ous << endl << "segments\n";
   printSegment(ous) << endl;
}
}
