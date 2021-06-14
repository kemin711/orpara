#include <sstream>

#include "interval.h"

namespace orpara {

/////// Interval Base class ///////////////////
Interval& Interval::operator=(const Interval& other) {
   if (this != &other) {
      b=other.b;
      e=other.e;
   }
   return *this;
}

bool Interval::less(const Interval &iv) const { 
   if (b < iv.b) return true;
   if (b > iv.b) return false;
   // now b and iv.b are the same, so compare e
   if (e < iv.e) return true;
   return false;
}

Interval::Interval(const string &str) {
   string::size_type i;
   i=str.find(",");
   if (i == string::npos) i=str.find("-");
   if (i == string::npos) i=str.find(";");
   if (i == string::npos) {
      cerr << "str is not accepted format for making an interval "
         << str << endl;
      exit(1);
   }
   b=atoi(str.substr(0,i).c_str());
   e=atoi(str.substr(i+1).c_str());
}

int Interval::overlap(const Interval &i) const
{
   if (isNull() && i.isNull()) return 0;
   if (isNull()) {
      return i.length();
   }
   if (i.isNull()) {
      return length();
   }
   if (e<i.b) return e-i.b+1; 
   else if (i.e<b) return i.e-b+1;
   else return min(e,i.e)-max(b, i.b) + 1;
}

int Interval::overlap(const int bb, const int ee) const {
   if (isNull()) return (ee-bb+1);
   if (e < bb) return e-bb+1;
   else if (ee < b)  return ee-b+1; 
   else return min(e,ee)-max(b,bb)+1;
}

int Interval::distance(const Interval& i) const {
   if (isNull() || i.isNull()) return 0;
   // [this interval] before [i interval]
   return max(b, i.b) - min(e,i.e);
}


/** base class merge does not invalidate input 
 * pointer.
 */
int Interval::merge(Interval *i, const int cut)
{
   if (isNull()) {
      b=i->b; e=i->e;
      return e-b+1;
   }
   int olp=overlap(*i);
   if (olp > cut) {
      b=min(b,i->b);
      e=max(e,i->e);
   }
   return olp;
}

int Interval::extend(const Interval &i) {
   if (isNull()) {
      b=i.b; e=i.e;
      return i.length();
   }
   int olp = overlap(i);
   b=min(b,i.b);
   e=max(e, i.e);
   return olp;
}


string Interval::toDelimitedString(const string& dl) const {
   ostringstream ost;
   ost << b << dl << e;
   return ost.str();
}

////// IntervalPile Class //////////

IntervalPile::~IntervalPile() {
   vector<Interval*>::iterator it;
   if (members.empty()) return;
   //cerr << "Calling IntervalPile destructor\n";
   //cerr << "member size: " << members.size() << endl;

   for (it=members.begin(); it != members.end(); it++) {
      delete *it;
   }
   members.clear();
}

/** does not actually destroy the elements in ip
 * when merging only playing with pointers.
 * Simply trnasfer ownership from ip to this object
 * for the underlying members.
 *
 * The members stays the same.
 */
int IntervalPile::merge(Interval *ip, const int cut) {
   int olp=Interval::merge(ip,cut);
   if (olp>cut)  {
      // transfer of ownership so that the members in ip
      // will not be destroyed
      IntervalPile* p=static_cast<IntervalPile*>(ip);
      members.insert(members.end(), p->members.begin(), p->members.end());
      //cerr << "after merge the member is nolonger valid\n";
      p->members.clear();
   }
   return olp;
}

void IntervalPile::add(Interval *i) {
   Interval::merge(i);
   members.push_back(i);
   i=0;
}

string IntervalPile::toDelimitedString(const string &dl) const {
   ostringstream ost;
   ost << Interval::toDelimitedString(dl) << dl 
      << size();
   return ost.str();
}

////////// HalfAlignPile Class //////
/** this method will only be called once for
 * each object
 * It will also compute the coverage.
 */
float HalfAlignPile::getIdentity() const {
   if (identity != -1.0) return identity;
   computeIntermediateResult();
   return identity;
}

void HalfAlignPile::computeIntermediateResult() const {
   vector<Interval*>::const_iterator i=members.begin();
   int sumiden=0;
   int sumalen=0;
   float sumcov=0;
   int sumlength=0;
   HalfAlignInterval *p;
   while (i != members.end()) {
      p = static_cast<HalfAlignInterval*>(*i);
      sumiden += p->getIdentical();
      sumalen += p->getAlnlen();
      sumcov += p->getCov();
      sumlength += p->length();
      ++i;
   }
   identity=(float)sumiden/sumalen;
   tcov=(float)sumcov/size();
   qcov=(float)sumlength/size()/length();
}

float HalfAlignPile::getQcov() const {
   if (qcov != -1.0) return qcov;
   computeIntermediateResult();
   /*
   vector<Interval*>::const_iterator i=members.begin();
   //float sumcov=0;
   int sumalen=0;
   HalfAlignInterval *p;
   while (i != members.end()) {
      p = static_cast<HalfAlignInterval*>(*i);
      //sumcov += p->getCov();
      sumalen += p->getAlnlen();
      ++i;
   }
   //cov=sumcov/size();
   cov=(float)sumalen/length();
   cov /= size();
   */
   return qcov;
}
float HalfAlignPile::getTcov() const {
   if (tcov == -1) {
      computeIntermediateResult();
   }
   return tcov;
}

string HalfAlignPile::toDelimitedString(const string &dl) const {
   ostringstream ost;
   ost << IntervalPile::toDelimitedString(dl)
      << dl << getIdentity() << dl 
      << getQcov() << dl << getTcov();
   return ost.str();
}

/////// IntervalChain Class ///////////////

IntervalChain::~IntervalChain() {
   //set<Interval*, lessInterval>::iterator i;
   //list<Interval*>::iterator i;
   //cerr << "Calling IntervalChain destructor: " << size() << " Nots\n";
   iterator i;
   for (i=begin(); i != end(); i++) 
      delete *i;
}

/** should not make copies of the Interval because
 * it could be any derived classes
 * Just adding pointers.
 */
void IntervalChain::add(Interval *iv) {
   //typedef set<Interval*, lessInterval>::iterator I;
   if (chain.empty()) {
      set(*iv);
      //chain.push_back(new Interval(*iv));
      chain.push_back(iv);
      return;
   }

   int cut=5;
   if (iv->getBegin() < getBegin()) 
      setBegin(iv->getBegin());
   if (iv->getEnd() > getEnd()) 
      setEnd(iv->getEnd());

   iterator i=begin();
   iterator ii, del;
   bool inserted=false;
   int olp;
   while (i != end()) {
      olp=(*i)->overlap(*iv);
      if (olp > cut) {
         //cerr << "Merging *i with iv\n";
         (*i)->merge(iv);
         inserted=true;
         delete iv;
         ii=i; ++ii;
         // consolidate possible 
         while (ii != chain.end() && !(*ii)->after(**i)) {
            if ((*i)->overlap(**ii) > cut) {
               //cerr << "merging *i and *ii\n";
               (*i)->merge(*ii); // transfer onwership if pile type
               del=ii;
               ++ii;
               delete *del;
               chain.erase(del);
            }
            else ++ii;
         }
         break;
      }
      else if (olp > 0 && olp <= cut &&  iv->before(**i)) {
         chain.insert(i, iv);
         inserted=true;
         delete iv;
         break;
      }
      else ++i;
   }
   if (!inserted) chain.push_back(iv); 
}

void IntervalChain::add(const int bb, const int ee) {
   //typedef set<Interval*, lessInterval>::iterator I;
   if (chain.empty()) {
      setBegin(bb);
      setEnd(ee);
      chain.push_back(new Interval(bb,ee));
      return;
   }
   if (bb < getBegin()) setBegin(bb);
   if (ee > getEnd()) setEnd(ee);

   //typedef list<Interval*>::iterator I;
   int cut=5;
   //I i=chain.begin();
   iterator i=chain.begin();
   //I ii, del;
   iterator ii, del;
   bool inserted=false;
   int olp;
   while (i != end()) {
      olp = (*i)->overlap(bb,ee);
      if (olp>cut) {
         (*i)->merge(new Interval(bb,ee));
         inserted=true;
         ii=i; ++ii;
         while (ii != chain.end() && !(*ii)->after(**i)) {
            if ((*i)->overlap(**ii) > cut) {
               (*i)->merge(*ii);
               del=ii;
               ++ii;
               delete *del;
               chain.erase(del);
            }
            else ++ii;
         }
         break;
      }
      else if (olp>0 && olp <= cut && bb < (*i)->getBegin()) {
         chain.insert(i, new Interval(bb, ee));
         inserted=true;
         break;
      }
      else ++i;
   }
   if (!inserted) chain.push_back(new Interval(bb,ee)); 
}

ostream& operator<<(ostream &ous, const IntervalChain &iv) 
{
   list<Interval*>::const_iterator i;
   for (i=iv.chain.begin(); i != iv.chain.end(); i++) {
      ous << (*i)->getBegin() << " - " << (*i)->getEnd() << endl;
   }
   return ous;
}

int IntervalChain::getInnerLength() const {
   list<Interval*>::const_iterator i=chain.begin();
   int suminner=0;
   while (i != chain.end()) {
      suminner += (*i)->length();
      ++i;
   }
   return suminner;
}

int IntervalChain::getOuterLength() const {
   return chain.back()->getEnd() - chain.front()->getBegin() + 1;
}

} // orpara namespace ends
