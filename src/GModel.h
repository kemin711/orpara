#ifndef GMODEL_H
#define GMODEL_H

#include <string>  
#include <vector> 
#include <list> 
#include <iostream>
#include <cstdlib>

/** class: GModel
 * Gene Model template class, no implementation.
 *
 * GModel.h This is an attemp to write a template function
 * on the ID of the genemodel.  The Id could be any type.
 *
 * I also implemented GenModel in the range
 * directory with different emphasis. 
 *
 * */
                
using namespace std; 

namespace orpara {
class ExondirectionError : public exception {
   public:
      ExondirectionError() : message("Different direction of exon inside the same gene") { }
      ExondirectionError(const string &msg) : message(msg) { }
      const char* what() const throw() { return message.c_str(); }
      ~ExondirectionError() throw() { }
   private:
      string message;
};
             
/** GModel
 * the class parameter is for the identifier of the model.
 * The ID could be either integer or string.
 *
 * A exon model of genes.
 * the exon is represented by a list of pair of integers.
 * This has much less functionality than the GModel 
 * derived from Noschain.  I need to rename one of them.
 *
 * This is an all header class, there is no corresponding cpp implementation
 * file.
 */       
template<typename T> class GModel;
template<typename T> 
ostream& operator<<(ostream& ous, const GModel<T> &mm);

template<class T> class GModel { 
   public: 
      GModel() : id(), exons(), intn(), intronsBuilt(false) { }
      GModel(const T &name) : id(name), exons(), intn(), intronsBuilt(false) { }
      GModel(GModel&& o) : id(std::move(o.id)), exons(std::move(o.exons)), intn(std::move(o.intn)), intronsBuilt(o.intronsBuilt) { }
      //void addExon(int bb, int ee) { exons.push_back(make_pair(bb,ee)); }
      
      /* require the exons are sorted from small to large
       */
      void addExon(int bb, int ee) throw(ExondirectionError);
      bool operator==(const GModel &mm) const;
      // return +, -, or ' ' if now known
      char direction() const;
      //pair<int,int> bound() const { return make_pair(exons[0].first, exons[exons.size()-1].second); }
      /* return the whole reange of the gene model
       */
      pair<int,int> bound() const { 
         return make_pair(exons.front().first, exons.back().second); }
      int size() const {
         if (direction() == '+')
            return exons.back().second - exons.front().first + 1; 
         else
            return exons.front().first - exons.back().second + 1; }
      /* one model contains another model
       * Both ends ousdide of mm, and the intron of 
       * mm is a subset of this object.
       */
      bool contain(const GModel &mm) const;
      /* only check for the introns, disregard the ends
       * This method is useful where some simple extension
       * of the partial models might have been wrong.
       */
      bool containIntron(const GModel &mm) const;
      /** 
       * checking only the begin and end of the model
       */
      bool endsContain(const GModel &mm) const;
      int numExons() const { return exons.size(); }
      //int start() const { return exons[0].first; }
      int start() const { return exons.front().first; }
      int begin() const { return exons.front().first; }
      //int end() const { return exons[exons.size()-1].second; }
      int end() const { return exons.back().second; }
      int finish() const { return exons.back().second; }
      // return all the intron boundaries
      const vector<int>& introns() const;
      // this is <T> strange <T> must be added after the function name.
      // output key, then exons
      friend ostream& operator<<<T>(ostream& os, const GModel<T>& mm);

      const T& getId() const { return id; }
      /* return all those introns longer than L
       */
      vector<pair<int,int> > intronLongerThan(const int L) const;
      int maxIntronLength() const;
      bool inside(pair<int,int> r) const {
         return min(start(), end()) > min(r.first, r.second) &&
            max(start(), end()) < max(r.first, r.second);
      }
      /** 
       * check whether two models overlap.
       * Requirement: they must be in the same direction.
       * Return: the length of the overlap.
       */
      int overlap(const GModel<T> &mm) const;
      /** 
       * Check overlap in opposite directions
       * Not implemented yet!
       */
      int overlapOpposite(const GModel<T> &mm) const { return 1; }
      /** 
       * This one should check the overlap regardless of direction
       * Not implemented yet.
       */
      int overlapBoth(const GModel<T> &mm) const { return 1; }
 
   protected: 
      T id; 
      //vector< pair<int,int> > exons;
      // sorted from 5'--to-->3' direction
      list< pair<int,int> > exons;
      mutable vector<int> intn;
      mutable bool intronsBuilt;
};  


template<class T>
void GModel<T>::addExon(int bb, int ee) throw(ExondirectionError) { 
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

template<class T>
bool GModel<T>::operator==(const GModel<T> &mm) const {
   if (exons.size() != mm.exons.size()) return false;
   list<pair<int,int> >::const_iterator i,j;
   for (i=exons.begin(), j=mm.exons.begin(); i != exons.end(); i++, j++) {
      if (i->first != j->first || i->second != j->second)
         return false;
   }
   return true;
}

template<class T>
char GModel<T>::direction() const {
   if (exons.empty()) return ' ';
   if (exons.front().first < exons.back().second)
      return '+';
   if (exons.front().first > exons.back().second)
      return '-';
   return ' ';
}

template<class T>
bool GModel<T>::contain(const GModel<T> &mm) const {
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

   int b=0, i;
   vector<int> Ithis = introns();
   vector<int> Ithat = mm.introns();
startSearch:
   while (b <= Ithis.size() - Ithat.size()) {
      //cout << "searching from " << b << endl;
      for (i=0; i<Ithat.size(); i++) {
         //cout << Ithat[i] << " x " << Ithis[b+i] << endl;
         if (Ithat[i] != Ithis[b+i]) {
            ++b;
            goto startSearch;
         }
      }
      /* debug code
      copy(Ithis.begin(), Ithis.end(), ostream_iterator<int>(cout, " "));
      cout << "\n contains intron: \n";
      copy(Ithat.begin(), Ithat.end(), ostream_iterator<int>(cout, " "));
      cout << endl;
      */
      return true;
   }
   /*
   cout << "\n+++++++++++++++++++++++\n";
   copy(Ithis.begin(), Ithis.end(), ostream_iterator<int>(cout, " "));
   cout << "\n does not contain intron: \n";
   copy(Ithat.begin(), Ithat.end(), ostream_iterator<int>(cout, " "));
   cout << endl;
   */
   return false;
}

// both gene must have at least one intron.
template<class T>
bool GModel<T>::containIntron(const GModel<T> &mm) const {
   if (numExons() < 2 || mm.numExons()<2
         || numExons() < mm.numExons()) return false;
   unsigned int b=0, i;
   vector<int> Ithis = introns();
   vector<int> Ithat = mm.introns();
startSearch:
   while (b <= Ithis.size() - Ithat.size()) {
      for (i=0; i<Ithat.size(); i++) {
         if (Ithat[i] != Ithis[b+i]) {
            ++b;
            goto startSearch;
         }
      }
      return true;
   }
   return false;
}

template<class T>
const vector<int>& GModel<T>::introns() const {
   if (intronsBuilt) return intn;
   if (!intn.empty()) intn.clear();
   if (numExons() == 1) return intn;

   //vector<int> tmp;
   list<pair<int,int> >::const_iterator i,j;
   i=exons.begin();
   j=--exons.end();
   if (direction() == '+') {
      intn.push_back(i->second+1);
      ++i;
      while (i != j) {
         intn.push_back(i->first-1);
         intn.push_back(i->second+1);
         ++i;
      }
      intn.push_back(j->first-1);
   }
   else if (direction() == '-') {
      intn.push_back(i->second-1);
      ++i;
      while (i != j) {
         intn.push_back(i->first+1);
         intn.push_back(i->second-1);
         ++i;
      }
      intn.push_back(j->first+1);
   }
   else {
      cerr << "direction wrong\n";
      exit(1);
   }
   intronsBuilt=true;
   return intn;
}

template<class T>
bool GModel<T>::endsContain(const GModel<T> &mm) const {
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

// out put the id first, then the exons
template<class T>
ostream& operator<<(ostream& ous, const GModel<T> &mm) {
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

template<class T>
vector<pair<int,int> > GModel<T>::intronLongerThan(const int L) const
{
   vector<pair<int,int> > tmp;
   vector<int> ins = introns();
   unsigned int i;
   if (direction() == '+') {
      for (i=0; i<ins.size(); i+=2) {
         if ((ins[i+1] - ins[i] + 1) > L) {
            tmp.push_back(make_pair(ins[i], ins[i+1]));
            //return true;
         }
      }
      //return false;
   }
   else if (direction() == '-') {
      for (i=0; i<ins.size(); i += 2) {
         if ( (ins[i]-ins[i+1]+1) > L) 
            tmp.push_back(make_pair(ins[i], ins[i+1]));
            //return true;
      }
      //return false;
   }
   //else return false;
   return tmp;
}

template<class T>
int GModel<T>::maxIntronLength() const {
   vector<int> ins = introns();
   int L=0;
   unsigned int i;
   if (direction() == '+') {
      for (i=0; i<ins.size(); i+=2) {
         if ((ins[i+1] - ins[i] + 1) > L) {
            L= ins[i+1] - ins[i] + 1;
         }
      }
   }
   else if (direction() == '-') {
      for (i=0; i<ins.size(); i += 2) {
         if ( (ins[i]-ins[i+1]+1) > L) 
            L= ins[i] - ins[i+1] + 1;
      }
   }
   return L;
}
template<class T>
int GModel<T>::overlap(const GModel<T> &mm) const
{
   if (direction() == '+') {
      if (mm.direction() == '+') {
         return min(end(),mm.end()) - max(begin(), mm.begin());
      }
      else if (mm.direction() == '-') {
         return 0;
      }
      else {
         cerr << "second model has no direction!\n";
         exit(1);
      }
   }
   else if (direction() == '-') {
      if (mm.direction() == '-') {
         return min(begin(), mm.begin()) - max(end(), mm.end());
      }
      else if (mm.direction() == '+') {
         return 0;
      }
      else {
         cerr << "second model has no direction!\n";
         exit(1);
      }
   }
   else {
      cerr << "first range has no direction()\n";
      exit(1);
   }
}
}
#endif
