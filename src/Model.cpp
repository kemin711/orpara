#include "Model.h"
#include <iterator>
#include <cstdlib>

namespace orpara {
//void Model::addExon(int bb, int ee) throw(ExondirectionError) { 
void Model::addExon(int bb, int ee) { 
   if (bb < ee) {
      if (numExons() > 0 && direction() == '-') {
         cerr << "must add exon to gene of the same direction\n"
            << " adding " << bb << "-->" << ee << " to "
            << *this << endl;
         throw ExondirectionError();
         //exit(1);
      }
      exons.push_back(make_pair(bb,ee)); 
   }
   else if (bb > ee) {
      if (numExons() > 0 && direction() == '+') {
         cerr << "must add exon to gene of the same direction\n"
            << " adding " << bb << "-->" << ee << " to "
            << *this << endl;
         throw ExondirectionError();
         //exit(1);
      }
      exons.push_front(make_pair(bb,ee));
   }
}

bool Model::operator==(const Model &mm) const {
   if (exons.size() != mm.exons.size()) return false;
   list<pair<int,int> >::const_iterator i,j;
   for (i=exons.begin(), j=mm.exons.begin(); i != exons.end(); i++, j++) {
      if (i->first != j->first || i->second != j->second)
         return false;
   }
   return true;
}

char Model::direction() const {
   if (exons.empty()) return ' ';
   if (exons.front().first < exons.front().second)
      return '+';
   if (exons.front().first > exons.front().second)
      return '-';
   return ' ';
}

bool Model::contain(const Model &mm) const {
   if (mm.numExons() > numExons() 
         || direction() != mm.direction()
         || !endsContain(mm)) 
      return false;
   // if the shorter one is single exon, 
   // needs special treatment
   if (mm.numExons()==1) {
      // both are single exon gene
      if (numExons() == 1)
         return true;
      list<pair<int,int> >::const_iterator it;
      if (direction() == '+') {
         for (it=exons.begin(); it != exons.end(); it++)
         {
            if (it->first <= mm.start() &&
                  it->second >= mm.end())
               return true;
         }
      }
      else if (direction() == '-') {
         for (it=exons.begin(); it != exons.end(); it++)
         {
            if (it->first >= mm.start() &&
                  it->second <= mm.end())
               return true;
         }
      }
      return false;
   }

   unsigned int b=0;
   vector<int> Ithis = introns();
   vector<int> Ithat = mm.introns();
startSearch:
   while (b <= Ithis.size() - Ithat.size()) {
      for (unsigned int i=0; i<Ithat.size(); i++) {
         if (Ithat[i] != Ithis[b+i]) {
            ++b;
            goto startSearch;
         }
      }
      return true;
   }
   return false;
}

// both gene must have at least one intron.
bool Model::containIntron(const Model &mm) const {
   if (numExons() < 2 || mm.numExons()<2
         || numExons() < mm.numExons()) return false;
   int b=0;
   vector<int> Ithis = introns();
   vector<int> Ithat = mm.introns();
startSearch:
   while (b <= signed(Ithis.size() - Ithat.size())) {
      for (unsigned int i=0; i<Ithat.size(); i++) {
         if (Ithat[i] != Ithis[b+i]) {
            ++b;
            goto startSearch;
         }
      }
      return true;
   }
   return false;
}

vector<int> Model::introns() const {
   vector<int> tmp;
   list<pair<int,int> >::const_iterator i,j;
   if (numExons() == 1) 
      return tmp;
   i=exons.begin();
   j=--exons.end();
   if (direction() == '+') {
      tmp.push_back(i->second+1);
      ++i;
      while (i != j) {
         tmp.push_back(i->first-1);
         tmp.push_back(i->second+1);
         ++i;
      }
      tmp.push_back(j->first-1);
   }
   else if (direction() == '-') {
      tmp.push_back(i->second-1);
      ++i;
      while (i != j) {
         tmp.push_back(i->first+1);
         tmp.push_back(i->second-1);
         ++i;
      }
      tmp.push_back(j->first+1);
   }
   else {
      cerr << "direction wrong\n";
      exit(1);
   }
   return tmp;
}

bool Model::endsContain(const Model &mm) const {
   if (direction() != mm.direction()) return false;

   if (direction() == '+') {
      if (start() <= mm.start() && end() >= mm.end())
         return true;
      return false;
   }
   else if (direction() == '-') {
      if (start() >= mm.start() && end() <= mm.end())
         return true;
      return false;
   }
   else return false;
}

ostream& operator<<(ostream& ous, const Model &mm) {
   ous << mm.id << endl;
   list<pair<int,int> >::const_iterator i;
   i=mm.exons.begin();
   ous << i->first << '-' << i->second;
   ++i;
   while (i != mm.exons.end()) {
      ous << "; " << i->first << '-' << i->second;
      ++i;
   }
   return ous;
}
}
