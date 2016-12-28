#include "FPGModel.h"
#include <cmath> // to solve the abs problem

namespace orpara {
ostream& operator<<(ostream &ous, const FPGModel &mm) {
   const char dl[]="\t";
   ous << mm.id << endl;
   ous << mm.qlen << dl << mm.tlen << dl << mm.qbegin << dl
      << mm.qend << dl << mm.tbegin << dl << mm.tend << dl
      << mm.numexon << dl
      << mm.sumscore << dl << mm.avgiden << dl << mm.qcov << dl 
      << mm.sumoverlap << endl;
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

// should return the number of covered
bool FPGModel::cover(const FPGModel &mm) const {
   if (numExons() < mm.numExons() || ! endsContain(mm)
         || getNormSumscore() < mm.getNormSumscore())
      return false;
   unsigned int covercnt=0, overlapcnt=0;
   if (direction() == '+' && mm.direction() == '+') {
      list<pair<int,int> >::const_iterator i,j;
      i=exons.begin();
      j=mm.exons.begin();
      while (i != exons.end() && j != mm.exons.end()) {
         while (i->second < j->first) ++i;
         if (i == exons.end()) break;
         else if (i->first <= j->first && i->second >= j->second) {
            ++covercnt; 
            ++i; ++j;
         }
         else if ((i->second >= j->first && i->second <= j->second)
               || (j->second >= i->first && j->second <= i->second))
         {
            //overlapping situation, and second cover first
             ++overlapcnt; ++i; ++j;
         }
         while (j != mm.exons.end() && j->second < i->first) ++j;
      }
   }
   else if (direction() == '-' && mm.direction() == '-') {
      // this is the mirror image of the code in the above section
      list<pair<int,int> >::const_reverse_iterator i,j;
      i=exons.rbegin();
      j=mm.exons.rbegin();
      while (i != exons.rend() && j != mm.exons.rend()) {
         while (i != exons.rend() && i->second > j->first) ++i;
         if (i == exons.rend()) break;
         else if (i->first >= j->first && i->second <= j->second) {
            ++covercnt; 
            ++i; ++j;
         }
         else if ((i->second <= j->first && i->second >= j->second)
               || (j->second <= i->first && j->second >= i->second))
         {
            // if second contain the first, then covered under
            // this overlapping situation.
             ++overlapcnt; ++i; ++j;
         }
         // sednod is before the first
         while (j != mm.exons.rend() && j->second > i->first)  ++j;
      }
   }
   else {
      cerr << "Error condition inside cover function\n";
      exit(1);
   }
   // make some judgment.
   if (covercnt == mm.exons.size()) {
      //cout << "All of mm.exon2 is covered\n";
      return true;
   }
   else {
      //cout << covercnt << " of " << mm.exons.size() << " covered\n"
       //  << overlapcnt << " overlap\n"; 
      return false;
   }
}

bool FPGModel::similar(const FPGModel &mm, float fr) const {
   if (direction() != mm.direction()
         || overlap(mm) < 50 || numExons() != mm.numExons()) 
      return false;
   double score1=getNormSumscore();
   double score2=mm.getNormSumscore();
   double iden1=getAvgiden();
   double iden2=mm.getAvgiden();
   /*
   if ( fabs(static_cast<double>( (score1-score2)/(0.5*(score1+score2)) )) < fr
         && fabs(static_cast<double>( (iden1-iden2)/(0.5*(iden1+iden2)) )) < fr 
         && fabs(static_cast<double>( (begin()-mm.begin())/(0.5*(begin()+mm.begin())))) < fr
         && fabs(static_cast<double>( (end()-mm.end())/(0.5*(end()+mm.end())))) < fr
         */
   if (  abs((score1-score2)/(0.5*(score1+score2)))  < fr
         && abs( (iden1-iden2)/(0.5*(iden1+iden2)) ) < fr 
         && abs( (begin()-mm.begin())/(0.5*(begin()+mm.begin()))) < fr
         && abs( (end()-mm.end())/(0.5*(end()+mm.end()))) < fr
       )
   {
      return true;
   }
   return false;
}
}
