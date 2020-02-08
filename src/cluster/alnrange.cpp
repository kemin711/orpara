#include "alnrange.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>


#ifdef HAVE_PG
// the query result must have fields in the 
// exact same order as 
//prt varchar(24), len integer, mb integer, me integer, score float, indetity float, ngidentity float, cov float);
//
alnrange::alnrange(PgDatabase* db, int i) {
	//prt=pgdb->GetValue(i,0);
	//len=atoi(pgdb->GetValue(i,1));
	begin=atoi(db->GetValue(i,2));
	end=atoi(db->GetValue(i,3));
	score=atof(db->GetValue(i,4));
	//identity=atof(pgdb->GetValue(i,5)); // not used yet
	ngidentity=atof(db->GetValue(i,6));
	cov=atof(db->GetValue(i,7));
}
#endif

float alnrange::ovlpfraction = 0.1;
int alnrange::ovlpcut = 5;

alnrange::alnrange(const alnrange &r) {
	begin=r.begin;
	end=r.end;
	score=r.score;
	ngidentity=r.ngidentity;
	cov=r.cov;
}

alnrange& alnrange::operator=(const alnrange& r) {
	if (this != &r) {
		begin=r.begin;
		end=r.end;
		score=r.score;
		ngidentity=r.ngidentity;
		cov=r.cov;
	}
	return *this;
}

bool alnrange::overlap(const alnrange &r, const int margin) {
	if (end < r.getBegin()+margin) return false;
	if (r.getEnd() < begin+margin) return false;
	return true;
}
bool alnrange::overlap(const alnrange &r) {
    int minrlen = end-begin+1 < r.end-r.begin+1 ? end-begin+1 : r.end-r.begin+1;
	if (end < r.getBegin()) return false;
	if (r.getEnd() < begin) return false;
    // now check ovlap situation
    int minend = end<r.end? end : r.end;
    int maxbegin = begin>r.begin? begin : r.begin;
    int overlap = minend - maxbegin;
    if (overlap < ovlpcut && overlap < minrlen*ovlpfraction) return false;
	return true;
}

string alnrange::asTabedString() const {
	ostringstream ost;
	ost << begin << '\t' << end << '\t'
		<< score << '\t' << ngidentity << '\t'
		<< cov;
	return ost.str();
}

string alnrange::asDelimitedString(const char sep[]) const {
	ostringstream ost;
	ost << begin << sep  << end << sep 
		<< score << sep << ngidentity << sep
		<< cov;
	return ost.str();
}
string alnrange::fields() {
	return "begin\tend\tscore\tngidentity\tcov";
}

string alnrange::essentialInfo() const {
	ostringstream ost;
	ost << begin << '\t' << end << '\t' << ngidentity;
	return ost.str();
}

///////////// rangePaire //////////////////
float rangePair::ngdiff_cut=0.25;
int rangePair::distance_cut=2000; // 50 kb apart
// for chlre2, avg=370, std=704, so 2000 is large enough
// min=1, max=34,245, this max intron size is not normal
// the mean+2*std intron size, this include 96% of the cases

rangePair::rangePair(const rangePair &r) : alnrange(r) {
	tprtid=r.tprtid;
	tcov=r.tcov;
	tmodelid=r.tmodelid;
	tgenomicid=r.tgenomicid;
	tstrand=r.tstrand;
	tstart=r.tstart;
	tend=r.tend;
}

string rangePair::asTabedString() const {
	ostringstream ost;
	ost << alnrange::asTabedString()
		<< '\t' << tprtid << '\t' << tmodelid << '\t'
		//<< tgenomicid << '\t' << tstrand << '\t'
		<< tstart << '\t' << tend << '\t'
		<< tcov;
	return ost.str();
}

// string variable is single-quoted
string rangePair::asDelimitedString(const char sep[]) const {
	string quote="";
	if (!strcmp(sep,"\t")) quote="";
	else if (!strcmp(sep,",")) quote="'";

	ostringstream ost;
	ost << alnrange::asDelimitedString()
		<< sep << quote << tprtid << quote << sep << tmodelid << sep
		//<< tgenomicid << '\t' << tstrand << '\t'
		<< tstart << sep << tend << sep
		<< tcov;
	return ost.str();
}

string rangePair::fields() {
	//return range::fields() + "\ttprtid\ttmodelid\ttgenomicid\ttstrand\ttstart\ttend\ttcov";
	return alnrange::fields() + "\ttprtid\ttmodelid\ttstart\ttend\ttcov";
}

string rangePair::essentialInfo() const {
	ostringstream ost;
	ost << alnrange::asTabedString()
		<< '\t' << tprtid << '\t' << tmodelid 
		<< '\t' << tstart << '\t' << tend;
	return ost.str();
}

/** this function is not very useful,
 * further more itoa is not implemented in 
 * this stdlib.h
list<string> rangePair::fieldsAsList() const {
	list<string> tmp;
	char buff[12];
	tmp.push_back(itoa(begin,buff));
	tmp.push_back(itoa(end,buff));
	tmp.push_back(gcvt(score,10,buff));
	tmp.push_back(gcvt(cov,3,buff));
	tmp.push_back(tprtid);
	tmp.push_back(itoa(tmodelid,10,buff));
	tmp.push_back(tgenomicid);
	tmp.push_back(tstrand);
	tmp.push_back(itoa(tstart,10,buff));
	tmp.push_back(itoa(tend,10,buff));
	tmp.push_back(gcvt(tcov,3,buff));
	return tmp;
}
*/

ostream& operator<<(ostream &ous, const alnrange& r) {
	ous << r.asTabedString();
	return ous;
	/*
	ous << r.begin << '\t' << r.end << '\t'
		<< r.score << '\t' << r.ngidentity << '\t'
		<< r.cov;
		*/
}

bool rangePair::sameGene(const rangePair &r) const {
	if (tgenomicid != r.tgenomicid || tstrand != r.tstrand) return false;
	if (abs(ngidentity - r.ngidentity) > ngdiff_cut) return false;
	if (tstrand == '+') {
		//cout << "on strand (+)\n";
		if (tend < r.tstart) {
			//cout << "correct order Left->Right\n";
			if (r.tstart-tend < distance_cut) {
				//cout << "within distance cut\n";
				return true;
			}
			else {
				//cout << "too far apart\n";
				return false;
			}
		}
		else {
			//cout << "incorrect order Right->Left\n";
			return false;
		}
	}
	else if (tstrand == '-') {
		//cout << "on strand (-)\n";
		//if (r.tend < tstart && tstart-r.tend < distance_cut) {
		if (r.tend < tstart) {
			//cout << "correct order <-R- <-L-\n";
			if (tstart-r.tend < distance_cut) {
				//cout << " with distance cut\n";
				return true;
			}
			else {
				//cout << "too far away\n";
				return false;
			}
		}
		else {
			//cout << "incorrect order on chrosome\n";
			return false;
		}
	}
	else {
		cerr << "strand must be + or - not " << tstrand << endl;
		exit(1);
	}
	return true;
}

//////////////////////// avgrange class ///////////////

// constructor, initialize the avgrange object
avgrange::avgrange(const alnrange &r) 
	: n(1), members(), covs(), sorted(false)
{
	sumcov = r.getCov();
	begin=r.getBegin();
	end=r.getEnd();
	sumbegin=r.getBegin();
	sumend=r.getEnd();
	sumscore=r.getScore()*r.getCov();
	sumng = r.getNg()*r.getCov();
	covs.push_back(r.getCov());  // this should be eliminated
	members.push_back(&r);
	//old methods, old code depends on this
	//method, do not add anything new
}

// rp must point to an object created on the heap
// rp=new range() or rp=new rangPair()
avgrange::avgrange(const alnrange *rp) 
	: n(1), members(), covs(), sorted(false)
{
	sumcov = rp->getCov();
	begin=rp->getBegin();
	end=rp->getEnd();
	sumbegin=rp->getBegin();
	sumend=rp->getEnd();
	sumscore=rp->getScore()*rp->getCov();
	sumng = rp->getNg()*rp->getCov();
	covs.push_back(rp->getCov());  // this should be eliminated
	//members.push_back(dynamic_cast<range*>(rp));
	members.push_back(rp);
}

avgrange::~avgrange() {
	vector<const alnrange*>::iterator it;
	for (it = members.begin(); it != members.end(); it++) {
		delete *it;
	}
	//members.clear();
}
// simple function, use the range class method
// for better control.
bool avgrange::overlap(const alnrange &r, const int margin) {
	if (end < r.getBegin()+margin) return false;
	if (r.getEnd() < begin+margin) return false;
	return true;
}

// having trouble with the polymorphism
// of range object
// the range must be a pointer to an object
// constructed by the outside world with the new operator
void avgrange::merge(const alnrange *r) {
	// do not merge the same object twice!
	for (size_t i=0; i<members.size(); i++) {
		if (r == members[i]) return;
	}
	
	++n;
	/* for reference format
	covs.push_back(r.getCov()); // more work
	sumbegin += r.getBegin();
	sumend += r.getEnd();
	sumscore += r.getScore()*r.getCov();
	sumng += r.getNg()*r.getCov();
	sumcov += r.getCov();
	begin = r.getBegin()<begin ? r.getBegin() : begin;
	end = r.getEnd()>end ? r.getEnd() : end;
	// how do we know the type of r?
	*/
	covs.push_back(r->getCov());
	sumbegin += r->getBegin();
	sumend += r->getEnd();
	sumscore += r->getScore()*r->getCov();
	sumng += r->getNg()*r->getCov();
	sumcov += r->getCov();
	begin = r->getBegin()<begin ? r->getBegin() : begin;
	end = r->getEnd()>end ? r->getEnd() : end;
	members.push_back(r);
}
void avgrange::merge(const avgrange *ar) {
	for (int i=0; i<ar->getCount(); i++) {
		//merge(ar->getMembers()[i]);
		merge(ar->members[i]);
	}
}

ostream& operator<<(ostream &ous, const avgrange &ar) {
	ous << "Range: " << ar.begin << " " << ar.end << endl
		<< "avg range: " << floor(ar.getBegin()) << " " << ceil(ar.getEnd())
		<< " N=" << ar.n << endl
		<< "avg: score=" << setprecision(0) << ar.getScore() 
		<< " ngidentity=" << setprecision(3) << ar.getNg() 
		<< " cov=" << setprecision(3) << ar.getCov()
		<< " median cov=" << ar.getMedianCov();
   return ous;
}

// write a table format
ostream& avgrange::writeTable(ostream &ous) const {
	//ous.unsetf(ios_base::scientific);
	//ous.setf(ios_base::fmtflags(0), ios_base::floatfield);
	ous << begin << "\t" << end << "\t";
	ous.setf(ios_base::fixed);
	ous << noshowpoint << setprecision(0)
		<< floor(getBegin()) << "\t" << ceil(getEnd())
		<< "\t" << n 
		<< "\t" << floor(getScore() + 0.5)
		<< "\t" << showpoint << setprecision(3) << getNg() 
		<< "\t" << setprecision(3) << getCov()
		<< "\t" << getMedianCov();
    return ous;
}
// for chimera table row output
string avgrange::asDelimitedString(const char sep[]) const {
	//ous.unsetf(ios_base::scientific);
	//ous.setf(ios_base::fmtflags(0), ios_base::floatfield);
	ostringstream ous;
	ous << begin << sep  << end << sep;
	ous.setf(ios_base::fixed);
	ous << noshowpoint << setprecision(0)
		<< floor(getBegin()) << sep << ceil(getEnd())
		<< sep << n 
		<< sep << floor(getScore() + 0.5)
		<< sep << showpoint << setprecision(3) << getNg() 
		<< sep << setprecision(3) << getCov()
		<< sep << getMedianCov();
    return ous.str();
}


ostream& avgrange::sqlinfo(ostream &ous) {
	ous.unsetf(ios::scientific);
	ous << "-- Range " << begin << '-' << end << " avg range " 
		<< floor(getBegin()) << '-' << ceil(getEnd()) << " N=" << n << endl
		<< "-- avg: score=" << setprecision(0) << getScore() 
		<< " ngidentity=" << setprecision(3) << getNg() 
		<< " cov=" << setprecision(3) << getCov()
		<< " median cov=" << getMedianCov();
	return ous;
}

double avgrange::getMedianCov() const { 
	if (!sorted) {
		sort(covs.begin(), covs.end());
		sorted = true;
	}
	if (n%2 == 0) {
		return 0.5*(covs[n/2-1] + covs[n/2]);
	}
	else return covs[n/2]; 
}

/** 1. ngidentity should be similar
 * comparison is directional, left to right
 */
string avgrange::checkSplit_debug(const avgrange &r) const
{
	string buff;
	//list<const range*>::const_iterator i,j;
	vector<const alnrange*>::const_iterator i,j;

	for (i=members.begin(); i != members.end(); i++) {
		for (j=r.members.begin(); j != r.members.end(); j++) {
			if ( static_cast<const rangePair*>(*i)->sameGene(*static_cast<const rangePair*>(*j)) ) 
			{
				buff += ((*i)->essentialInfo() 
							+ '\t' 
							+ (*j)->essentialInfo() 
							+ "\n");
			}
		}
	}
	return buff;
}

// this one is for real out put, more flexible
// check the genomic location information and difference in
// sequence identity
list<string> avgrange::checkSplit(const avgrange &r, const char sep[]) const
{
	//list<string> buff;
	list<string> result;
	//list<const range*>::const_iterator i,j;
	vector<const alnrange*>::const_iterator i,j;

	for (i=members.begin(); i != members.end(); i++) {
		for (j=r.members.begin(); j != r.members.end(); j++) {
			if ( static_cast<const rangePair*>(*i)->sameGene(*static_cast<const rangePair*>(*j)) ) 
			{
				result.push_back( "'" + static_cast<const rangePair*>(*i)->genomicInfo() + "'" 
						+ sep 
						+ (*i)->asDelimitedString(sep) + sep
						+ (*j)->asDelimitedString(sep));
			}
		}
	}
	return result;
}

SplitResult avgrange::testSplit(const avgrange &r, const char sep[]) const
{
	//list<SplitResult> result;
	SplitResult res;
	vector<const alnrange*>::const_iterator i,j;

	for (i=members.begin(); i != members.end(); i++) {
		for (j=r.members.begin(); j != r.members.end(); j++) {
			if ( static_cast<const rangePair*>(*i)->sameGene(*static_cast<const rangePair*>(*j)) ) 
			{
				res.add(static_cast<const rangePair*>(*i), static_cast<const rangePair*>(*j));
			}
		}
	}
	return res;
}

list<string> SplitResult::outputRow(const char sep[]) const {
	list<string> res;
	string quote="";
	if (!strcmp(sep,"\t")) quote="";
	else if (!strcmp(sep,",")) quote="'";

	list< pair<rangePair,rangePair> >::const_iterator it;
	for (it=joins.begin(); it != joins.end(); it++) {
		ostringstream ous;
		ous << quote << guide << quote << sep << guideLen 
			<< sep << quote << it->first.getGenomic() << quote 
			<< sep << quote << it->first.getStrand() << quote
			<< sep << it->first.asDelimitedString(sep) 
			<< sep << it->second.asDelimitedString(sep);
		res.push_back(ous.str());
	}
	return res;
}
