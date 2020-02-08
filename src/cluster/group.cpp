#include "group.h"
#include <stdlib.h>
#include <iomanip>
#include <math.h>

/* The feature order is hard-coded into this program.  Not a good
 * design. Input data should be in this order. */
const char* divstat::features[FEAT] = { 
   "qcov", "tcov", "score", "identity",
		"ngidentity", "similarity", "matchlen", "nogaplen" };

/* initialize the vector to zero */
divstat::divstat() {
	id[0] = '\0';
	for (int i=0; i<FEAT; i++) feat[i] = 0;
}

divstat::divstat(istream &in) {
	in >> id;
	for (int i=0; i<FEAT; i++) in >> feat[i];
}
 
void divstat::read(istream &in) {
	in >> id;
	for (int i=0; i<FEAT; i++) in >> feat[i];
}
/* must select all features as defined in features array.  The order
 * of the selection is not essential.  This is a little slower than
 * if the order of the columns are known. */
//divstat::divstat(PgDatabase &pgdb, int row) {
divstat::divstat(result &qres, int row) {
	if (qres.size()<FEAT) {
		cerr << "Not enought columns selected from the database\n";
		for (int i=0; i<FEAT; i++) {
			cerr << features[i] << " ";
		}
		cerr << endl << " are required\n";
	}
	for (int i=0; i<FEAT; i++) {
		//feat[i]=atof(pgdb.GetValue(row, features[i]));
		//feat[i]=atof(qres[row][features[i]]);
		qres[row][features[i]].to(feat[i]);
	}
}
	
divstat::divstat(const divstat &d) {
	strcpy(id, d.id);
	for (int i=0; i<FEAT; i++) feat[i] = d.feat[i];
}

/* nice output for human, may not work for database input */
ostream& operator<<(ostream &o, const divstat &ds) {
	o << setw(5) << ds.id << '\t';
	o.setf(ios::showpoint);
	o.setf(ios::fixed, ios::floatfield);
	o.precision(3);
	o << '\t' << ds.feat[0] << '\t' << ds.feat[1] << '\t';
	o.unsetf(ios::showpoint);
	o << setprecision(0) << ds.feat[2] << '\t'; 
	o.setf(ios::showpoint);
	o << setprecision(3) << ds.feat[3] << '\t'
		<< ds.feat[4] << '\t' << ds.feat[5] << '\t';
	o.unsetf(ios::showpoint);
	o << setprecision(0) << ds.feat[6] << '\t' << ds.feat[7];
	return o;
}
ostream& divstat::dump(ostream &ou) {
	ou << id;
	for (int i=0; i<FEAT; i++) ou << '\t' << feat[i];
	return ou;
}

void divstat::scale(double f) {
	for (int i=0; i<FEAT; i++) feat[i] *= f;
}
divstat& divstat::operator=(const divstat &d) {
	if (this != &d) {
		strcpy(id, d.id);
		for (int i=0; i<FEAT; i++) feat[i] = d.feat[i];
	}
	return *this;
}
bool divstat::operator>(const divstat &ds) const { 
	//return (iden>ds.iden && mlen>ds.mlen && 
	//	qcov>ds.qcov && tcov>ds.tcov && sc>ds.sc); 
	//	only count score and idenity
	//return (iden>ds.iden && qcov>ds.qcov && sc>ds.sc); 
	//return (feat[0]>ds.feat[0] && feat[2]>ds.feat[2] && feat[4]>ds.feat[4]); 
	return (getScore() > ds.getScore());
}
bool divstat::operator<(const divstat &ds) const { 
//	return (iden<ds.iden && mlen<ds.mlen && 
//		qcov<ds.qcov && tcov<ds.tcov && sc<ds.sc); 
	//return (iden<ds.iden && qcov<ds.qcov && sc<ds.sc); 
	//return (feat[0]<ds.feat[0] && feat[2]<ds.feat[2] && feat[4]<ds.feat[4]); 
	return (getScore()<ds.getScore());
}
bool divstat::smallerThan(const divstat &d, double k)  const {
	if (!&d) return false;  // this can save a lot of testing
	divstat tmp = d;
	tmp.scale(k);
	return (*this < tmp);
}
bool divstat::smallerThan(const divstat *const d, const double k) const {
	if (!d) return false;
	divstat tmp = *d;
	tmp.scale(k);
	return (*this < tmp);
}
/***********   end of divstat  ****************************************/
////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//////////////         group class           /////////////////////////
/////////////////////////////////////////////////////////////////////

/* initialize the object to impossible values 
 * sapce for divm vector not allocated yet. 
 * */
group::group() : query(-1), divCnt(0), anchor(-1) {
	if (divisions.empty()) {
		cerr << "divisions not initialized yet\n";
		exit(1);
	}
	divm.resize(divisions.size());
	int i;
	// no memory allocated for dimv vector. It is the job of next.
	for (i=0; i<divisions.size(); i++) divm[i] = 0;
	for (i=0; i<FEAT; i++) {
		stat[i][0] = 0;
		stat[i][1] = -1;
	};
}

/* copy constructor */
group::group(const group &g) : divm(divisions.size()){
	query = g.query;
	for (int i=0; i<divisions.size(); i++) {
		if (g.divm[i]) {
			divm[i] = new divstat;
			*divm[i] = *g.divm[i];
		}
		else divm[i] = 0;
	}
	divCnt = g.divCnt;
	anchor = g.anchor;
	for (int i=0; i<FEAT; i++) {
		stat[i][0] = g.stat[i][0];
		stat[i][1] = g.stat[i][1];
	}
}
/* the initialization, or input constructor 
 * it is used in the training program.  */
group::group(int &q, istream &in) throw(inputend) {
	divm.resize(divisions.size());
	for (int i=0; i<divisions.size(); i++) {
		divm[i] = 0;
	}
	divCnt = 0;
	query = q;
	do {
		divCnt++;
		divstat* tmp = new divstat(in);
		divm[getdividx(tmp->getid())] = tmp;
		in >> q;
	} while (!in.eof() && q == query);
	// set up anchor an anchor should always be found
	anchor = -1;
	nextAnchor();
	if (anchor == -1) { cerr << "anchor not found\n"; exit(1); } // debug
	dostat();  // initialize the stat member
	if (in.eof()) {
		cout << "input ended\n";  // debug
		throw inputend();
	}
}

// behaves like a constructor, default constructor already
// allocated space, this function needs to fill up the members
bool group::next(int &q, istream &in) {
	// done in default constructor divm.resize(divisions.size());
	for (int i=0; i<divisions.size(); i++) {  // clean up routine
		delete divm[i];
		divm[i] = 0;
	}
	divCnt = 0;
	query = q;
	do {
		divCnt++;
		divstat* tmp = new divstat(in);
		divm[getdividx(tmp->getid())] = tmp;
		in >> q;
	} while (!in.eof() && q == query);
	// set up anchor. an anchor should always be found
	anchor = -1;
	nextAnchor(); // seearch from (old anchor + 1)
	if (anchor == -1) { cerr << "anchor not found\n"; exit(1); } // debug
	dostat();  // initialize the stat member
	if (in.eof()) {
		cout << "input ended\n";  // debug
		//throw inputend();
		return false;
	}
	return true; // more data following
}
group& group::operator=(const group &g) {
	if (this != &g) {
		query = g.query;
		divCnt = g.divCnt;
		anchor = g.anchor;
		for (int i=0; i<divisions.size(); i++) {
			if (!divm[i]) { // null for this divstat
				if (g.divm[i]) {  // input not null, allocate space and copy
					divm[i] = new divstat;
					*divm[i] = *g.divm[i];
				}
			}
			else { // this not null
				if (!g.divm[i]) {  // input null, needs deallocate memory
					delete divm[i];
					divm[i] = 0;
				}
				else *divm[i] = *g.divm[i];
			}
		}
		/* this part is not needed, wasted */
		for (int i=0; i<FEAT; i++) {
			stat[i][0] = g.stat[i][0];
			stat[i][1] = g.stat[i][1];
		}
	}
	return *this;
}

// set the internal anchor member value
int group::nextAnchor() {
	int i = anchorIndex(anchor) + 1;
	while (i < divisions.size() && !divm[ searchOrder[i] ] ) i++;
	if (i == divisions.size()) return -1;
	anchor = searchOrder[i];
	return anchor;
}

/* output raw data for human to read */
ostream& operator<<(ostream &o, const group &g) {
	o << "Query: " << g.query << endl;
	for (int i=0; i<g.divisions.size(); i++) {
		o << '\t';
		if (g.divm[i]) o << *(g.divm[i]) << endl;
		else o << "NULL" << endl;
	}
	return o;
}

void group::rmbranch(int i, ostream &ou) {
	ou << *divm[i] << "branch discarded from" << query << "\n"; //debug
	rmbranch(i);
}
void group::rmbranch(int i) {
	delete divm[i];
	divm[i] = 0;
	divCnt--;
	// if anchor is the one to be removed
	if (i == anchor) nextAnchor();

	if (anchor == -1) {   // debug stage
		cerr << "anchor set to -1\n";
		exit(1);
	}

	dostat();  // initialize the stat member
}

/* return true if branch i is revomed */
bool group::checkAndRemove(int i, int b, int e, ostream &ou) {
	for (int j=b; j<e; j++) {
		if (divm[j] && *divm[j] > *divm[i]) {
			rmbranch(i, ou);
			return true;
		}
	}
	return false;
}

pair<int, double> group::mintcov() const {
	int i, mini;
	double min = 100000;
	for (i=0; i<divisions.size(); i++) {
		if (divm[i] && divm[i]->getTcov() < min) {
			min = divm[i]->getTcov();
			mini = i;
		}
	}
	return make_pair(mini, min);
}

double group::getMinTcov(int &div) const {
	double tmp = 100000;
	for (int i=0; i<divisions.size(); i++) {
		if (divm[i] && divm[i]->getTcov() < tmp) {
			tmp = divm[i]->getTcov();
			div = i;
		}
	}
	if (tmp > 1) {
		cerr << "Inside getMinTcov: minimum of coverage larger than 1\n";
		exit(1);
	}
	return tmp;
}
pair<int, double> group::max(int f) const {
	int maxi, i;
	double max = 0;
	for (i=0; i<divisions.size(); i++) {
		if (!divm[i]) continue;
		if (divm[i]->feat[f] > max) {
			max = divm[i]->feat[f];
			maxi = i;
		}
	}
	//cout << "max " << divstat::features[f] << " is " << max 
	//	<< " div " << divisions[maxi] << " got the max\n";
	return make_pair(maxi, max);
}


// output the raw data
void group::dumpAsTable(ostream &ou) {
	for (int i=0; i<divisions.size(); i++) 
		if (divm[i]) ou << query << '\t' << *divm[i] << endl;
}
void group::dumpKey(ostream &ou) {
	for (int i=0; i<divisions.size(); i++) 
		if (divm[i]) ou << query << '\t' << divm[i]->id << endl;
}
ostream& group::dump(ostream &ou) {
	for (int i=0; i<divisions.size(); i++)
		if (divm[i]) {
			ou << query << '\t';
			divm[i]->dump(ou) << endl;
		}
	return ou;
}
void group::dumpStat(ostream &ou) {
	for (int i=0; i<FEAT; i++) 
		ou << stat[i][0] << "," << stat[i][1] << '\t';
}
/* statistics of each feature over all matching target divisions
 * if the sequence length is about 40000 (longest protein),
 * the formula sigma X^2 may produce overlow problem.
 * */
void group::dostat() {
	//double stat[5][2] = {0}; //for iden, mlen, qcov, tcov, sc
	//double sumRatio = 0;
	/* use a different formula for computing mean and standard deviation
	 * The recursive formula that has not overflow problem.
	 * M1 = X1 
	 * V1 = 0
	 * Mi+1 = Mi + (Xi+1 - Mi)/(i+1)
	 * Si+1 = (1 - 1/i)Si + (i+1)(Mi+1 - Mi)^2  // sample variance
	 * Vi+1 = j(Vj/(j+1) + (Xj+1 - Xj)^2)
	 * */

	int g=0;  // group index
	while (!divm[g]) g++;  // the first non empty target division
	int f;  // features
	for (f=0; f<FEAT; f++) {
		stat[f][0] = divm[g]->feat[f];    // mean
		stat[f][1] = 0;    // Variance
	}  // initialize
	int i;
	int cnt = 0;
	double Xi;
	for (i=++g; i<divisions.size(); i++) {
		if (divm[i]) {
			cnt++;
			for (f=0; f<FEAT; f++) {
				Xi = stat[f][0];  // tmp variable
				stat[f][0] = Xi + (divm[i]->feat[f] - Xi)/(cnt+1);
				// now stat[f][0] is Xi+1
				// Sample Variance Formula
				// stat[f][1] = (1-1/cnt)*stat[f][1] + (cnt+1)*pow((stat[f][0] - Xi), 2);
				stat[f][1] = cnt*(stat[f][1]/(cnt+1) + pow((stat[f][0] - Xi), 2));
			}
		}
	}
	for (i=0; i<FEAT; i++) {
		stat[i][1] = sqrt(stat[i][1]);
		//cout << stat[i][0] << ", " << sqrt(stat[i][1]) << '\t';
	}
}
void group::dumpAllSMratio(ostream &ou) {
	for (int i=0; i<FEAT; i++) {
		ou << stat[i][1]/stat[i][0] << " ";
	}
}

bool group::isConserved() const {
	if (stat[divstat::idxidentity()][0] > csvIdentity) {
		if (stat[divstat::idxidentity()][1] < csvStd ||
				stat[divstat::idxngidentity()][1] < csvStd) {
			if (getScoreSMratio() < 0.1) return true;
			else return false;
		}
		else if (highVarDivPresent() && 
				(stat[divstat::idxidentity()][1] < csvStdHighVar ||
				 stat[divstat::idxngidentity()][1] < csvStdHighVar)) {
			if (getScoreSMratio() < 0.1) return true;
			else return false;
		}
	}
	return false;
}

bool group::highVarDivPresent() const{
	if (divm[gconst::getdividx(highVarDiv[0])] || divm[gconst::getdividx(highVarDiv[1])])
		return true;
	return false;
}

///////////////////////////////////////////////////////////////////
////////////////  gdiagnosis class   //////////////////////////////

/*** initialize static members 
 * the Guid matrix has divisions.size() rows
 *                     FEAT columns with two vallues(avg, std)
 * Because divisions must have been initialized, we can use it safely.
 **************/

vector< vector< vector<double> > > gdiagnosis::guid=
		vector< vector< vector<double> > >(divisions.size());

/* We are reading from the output from the training class: trainer
 * */
void gdiagnosis::readGuid(istream &in) {
	// comments starts with #
	string dumy; // first two lines are comments
	char ch;
	in.get(ch);
	while (ch == '#') {     // starting a comment line
		getline(in, dumy);
		in.get(ch);
	}
	guid.resize(divisions.size());

	for (int i=0; i<divisions.size(); i++) {
		guid[i].resize(FEAT);
		in >> dumy;   // discarding div name, fixed order; known to program
		for (int j=0; j<FEAT; j++) {
			guid[i][j].resize(2);
			in >> guid[i][j][0] >> guid[i][j][1];
		}
		in >> dumy;  // discarding count
	}
}
///////////////// end of static member functions ///////////////

gdiagnosis::gdiagnosis() 
	: group(), zcut(1)
{
	getspace();
}

gdiagnosis::gdiagnosis(double cut) 
	: group(), zcut(cut)
{
	getspace();
}

/* not used any more
gdiagnosis::gdiagnosis(int &q, istream &in) throw(inputend) 
try  	: group(q, in) 
{ 
	getspace();
	calzval();
} 
catch (inputend &x) { throw x; } 
*/
bool gdiagnosis::next(int &q, istream &in) 
{ 
	bool reachedEnd = group::next(q,in);
	//getspace();  no needs to zero all entries
	// space already there in constructed objects
	calzval();
	return reachedEnd;
} 

/* copy constructor, reallocate space for zval and norm vector */
gdiagnosis::gdiagnosis(const gdiagnosis &gd) 
	: group(gd) 
{
	getspace();
	int i,j;
	for (i=0; i<divisions.size(); i++) {
		for (j=0; j<FEAT; j++) {
			zval[i][j] = gd.zval[i][j];
			norm[i][j] = gd.norm[i][j];
		}
		zval[i][FEAT] = gd.zval[i][FEAT];
	}
	for (i=0; i<FEAT; i++) {
		zval[divisions.size()][i] = gd.zval[divisions.size()][i];
		for (j=0; j<4; j++) hgrm[j][i] = gd.hgrm[j][i];
	}
	zval[divisions.size()][FEAT] = gd.zval[divisions.size()][FEAT];
}

gdiagnosis& gdiagnosis::operator=(const gdiagnosis& g) {
	if (this != &g) {
		this->group::operator=(g);
		int i,j;
		for (i=0; i<divisions.size(); i++) {
			for (j=0; j<FEAT; j++) {
				zval[i][j] = g.zval[i][j];
				norm[i][j] = g.norm[i][j];
			}
			zval[i][FEAT] = g.zval[i][FEAT]; // last column of zval
		}
		// copy buttom row, and the hgrm matrix
		for (i=0; i<FEAT; i++) {
			zval[divisions.size()][i] = g.zval[divisions.size()][i];
			for (j=0; j<4; j++) hgrm[j][i] = g.hgrm[j][i];
		}
		zval[divisions.size()][FEAT] = g.zval[divisions.size()][FEAT];
	}
	return *this;
}

/* used the anchor as standard, all values are 10 
 * missing divisions are assigned zero values*/
void gdiagnosis::calzval() {
	if (divCnt < 2) return;  // only do this if more than 1 div
	// zval, hgrm not changed
	int i, j;

	// buttom row set to zero; hgrm set to zero
	for (j=0; j<FEAT+1; j++) {
		zval[divisions.size()][j] = 0;
		for (i=0; i<4; i++) hgrm[i][j] = 0;
	}
	for (i=0; i<divisions.size(); i++) {
		if (divm[i]) { // exists
			zval[i][FEAT] = 0;  // last column set to zero
			if (i == anchor) {
				for (j=0; j<FEAT; j++) {
					norm[i][j] = guid[i][j][0]; // copy the value from guid
					zval[i][j] = 0;
				}
			}
			else {  // not an anchor, needs more computation
				for (j=0; j<FEAT; j++) {
					norm[i][j] = guid[anchor][j][0]*divm[i]->feat[j]/divm[anchor]->feat[j];
					zval[i][j] = (norm[i][j] - guid[i][j][0])/guid[i][j][1];

					// fill in the hgrm matrix
					if (zval[i][j] < -zcut) hgrm[0][j]++;
					else if (zval[i][j] > zcut) hgrm[2][j]++;
					else hgrm[1][j]++;
					hgrm[3][j]++;  // total count

					// buttom row accumulate
					zval[divisions.size()][j] += fabs(zval[i][j]);
					zval[i][FEAT] += zval[i][j]; // last column 
				}
			}
			zval[divisions.size()][FEAT] += zval[i][FEAT];
		}
		else for (j=0; j<FEAT; j++) {
			zval[i][j] = 0; // set to zero, not interfere with calculation
		}
	}
	
	// no need to average the buttom row any more
	for (i=0; i<FEAT; i++) { // average the buttom row
		zval[divisions.size()][i] = zval[divisions.size()][i]/(divCnt-1);
	}
}

/* If at least one division is marked as removed, return true
 * If the anchor is bad, this function will fix it.  Non-anchor
 * divisions will simply be marked as removed.
 * */
vector<int> gdiagnosis::rmlow(double lowercut, double avgzc) {
	int i, j;
	int iq, ii, is, ing, isc;  // index
	iq = divstat::idxqcov();
	ii = divstat::idxidentity();
	is = divstat::idxsimilarity();
	ing = divstat::idxngidentity();
	isc = divstat::idxscore();
	bool pqc, pid, psc;  // passed qc (qcov) id (identities), sc (score)
	vector<int> poordiv;

tryagain: poordiv.clear();   // retry after removing anchor

	for (i=0; i<divisions.size(); i++) {
		if (!divm[i]) continue;
		if (i == anchor) {
			if (avgzval() > avgzc) {  // anchor quality questionable
				if (anchorPartial()) { // anchor partial
					// needs to use the next anchor and recalculate  everything
					nextAnchor();
					if (anchor == -1) {  // debug
						cerr << "Running out of anchor " << query << endl;
						exit(1);
					}
				}
				else {  // truly remove, also needs to recalculate everything
					rmbranch(i);
					poordiv.push_back(i);
					if (divCnt < 2) break;
				}
				calzval();
				goto tryagain;
			}
		}
		else {
			pqc=true, pid=true, psc=true;
			if (zval[i][iq] < -zcut) pqc = false;
			if (zval[i][isc] < -zcut) psc = false;
			if (zval[i][ii] < -zcut && zval[i][is] < -zcut && 
					zval[i][ing] < -zcut) pid = false;

			if ( !(pqc && psc && zval[i][ii] > lowercut && zval[i][is] > lowercut && zval[i][ing] > lowercut) &&
					!(psc && pid && zval[i][iq] > lowercut) &&
					!(pid && pqc && zval[i][isc] > lowercut) ) 
			{
				rmbranch(i);
				poordiv.push_back(i);
				if (divCnt < 2) break;
				calzval();
			}
		}
	}

	return poordiv;
}

/////////////////  diagnostic functions ////////////////////////////////
bool gdiagnosis::qualitypass() { 
	int i, ngi, s;
	i =divstat::idxidentity();
	ngi =divstat::idxngidentity();
	s =divstat::idxsimilarity();

	return (hgrm[0][i] == 0 && hgrm[2][i] != hgrm[3][i]) ||
	         (hgrm[0][ngi] == 0  && hgrm[2][ngi] != hgrm[3][ngi]) || 
				(hgrm[0][s] == 0 && hgrm[2][s] != hgrm[3][s]); 
}

/* Average fo zval for all featues smaller than zcut */
bool gdiagnosis::passed(const double zcut) const {
	for (int i=0; i<FEAT; i++) {
		if (zval[divisions.size()][i] > zcut) return false;
	}
	return true;
}
bool gdiagnosis::passedAvg(const double zcut) const {
	//double sum=0;
	//for (int i=0; i<FEAT; i++) sum += zval[divisions.size()][i];
	return (getAvgZval() < zcut);
}
bool gdiagnosis::goodQuality(double zcut) const {
	if (divCnt == 1 || isConserved() || passedAvg(zcut)) return true;
	return false;
}
bool gdiagnosis::lastTest() const {
	if (zval[divisions.size()][divstat::idxidentity()] < 1 && zval[divisions.size()][divstat::idxscore()] < 1) return true;
	else if (zval[divisions.size()][divstat::idxscore()] < 1.5) return true;
	return false;
}
bool gdiagnosis::passedQcov(const double zcut) const {
	return zval[divisions.size()][divstat::idxqcov()] < zcut;
}
bool gdiagnosis::passedIdentity(const double zcut) const {
	return (zval[divisions.size()][divstat::idxidentity()] < zcut ||
		zval[divisions.size()][divstat::idxngidentity()] < zcut ||
		zval[divisions.size()][divstat::idxsimilarity()] < zcut);
}
bool gdiagnosis::passedLength(const double zcut) const {
	return (zval[divisions.size()][divstat::idxmatchlen()] < zcut ||
		zval[divisions.size()][divstat::idxnogaplen()] < zcut);
}

/* needs to get the philocerphy of this function */
bool gdiagnosis::coverageDefect(double zcut) const {
	/* identity and qcov close to middle */
	if (zval[divisions.size()][divstat::idxidentity()] < zcut && zval[divisions.size()][divstat::idxtcov()] < zcut && 
			getTcovSMratio() < 0.2) {
		return true;
	}
	else if (getIdenSMratio() < 0.07 && stat[divstat::idxidentity()][0] > 50 && getMlenSMratio() < 0.07 && getQcovSMratio() < 0.07 && stat[divstat::idxqcov()][0] > 0.5) {
		return true;
	}
	else if (getTcovSMratio() < 0.1 && stat[divstat::idxtcov()][0] > 0.5 && getIdenSMratio() < 0.1) {
		return true;
	}
	return false;
}

///////////////////  helper functions ///////////////////////////

/* try to remove branches that do not belong, then do a test
 * of quality, if good return true */
bool gdiagnosis::trimAndTest(ostream &ou, double zcut) {
	/* anchor too low */
	if (zval[divisions.size()][FEAT] > 0) {  // larger than anchor
		if (zval[divisions.size()][FEAT]/(divCnt-1) > 7.5) {
			rmbranch(anchor, ou);     // debug info
			calzval();
			if (goodQuality(zcut) || coverageDefect(0.8)) return true;
		}
	}
	else { // anchor OK, none anchor bad
		bool trimmed = false;
		for (int i=0; i<divisions.size(); i++) {
			if (divm[i] && zval[i][FEAT] < -7.5 ) {
				rmbranch(i, ou);
				trimmed = true;
			}
		}
		if (trimmed) {
			calzval();       // only needs this if this object has been
			if (goodQuality(zcut) || coverageDefect(0.8)) return true;
		}
	}
	return false;
}
/* this method is not universal 
 * based on hs->fr identity should be smaller than hsl->non-fish vertebrates
 * */
void gdiagnosis::trimByIden() {
	if (divm[6]) { // fr divisions exist
		for (int i=0; i<6; i++) {
			if (divm[i] && divm[i]->feat[0] < divm[6]->feat[0]) {
				rmbranch(i, cerr);
				calzval();
			}

		}
	}
	else if (divm[7]) { // vrt exist
		for (int i=0; i<6; i++) {
			if (divm[i] && divm[i]->feat[0] < divm[7]->feat[0]) {
				rmbranch(i, cerr);
				calzval();
			}
		}
	}
	else if (divm[0]) { // pri
		for (int i=1; i<8; i++) {
			if (divm[i] && divm[i]->feat[0] > divm[0]->feat[0]) {
				rmbranch(0, cerr);
				calzval();
				break;
			}
		}
	}
}


bool gdiagnosis::fixScoreAndTest(double szcut, gdiagnosis &d) const {
	/* first check whether it is fixable or not, then do the fixing
	 * and checking */
	//if (zval[8][0] > 1.49) return false;
	if (zval[divisions.size()][0] < 1.5) {  // which ones are we testing
		d = *this;
		for (int i=0; i<divisions.size(); i++) {
			if (divm[i] && divm[i]->feat[divstat::idxqcov()] < 1) {
				d.divm[i]->feat[divstat::idxscore()] /= d.divm[i]->feat[divstat::idxqcov()];
			}
		}
		d.calzval();
		if (d.zval[divisions.size()][divstat::idxscore()] < szcut) return true;
	}
	return false;
}

bool gdiagnosis::targetPartial(int i) const { // if only the division is not an anchor
	if (i == anchor) { 
		cerr << "I have not considered anchor yet\n";
		exit(1);
	}
	// try to correct score, mlen, and nogaplen based on qcov
	// if meets the basic test
	const double IR = 0.81;  // for div identity > IR*anchor_ideitity
	const double IZC = -0.5;  // zcut for identity, nogapidentity, and
									 // similarity consideration
	

	double qcov, anchor_qcov, mlen, anchor_mlen, nglen, anchor_nglen; 
	double score;
	//double tcov, anchor_tcov;
	qcov = divm[i]->getQcov();
	anchor_qcov = divm[anchor]->getQcov();
	mlen = divm[i]->getMatchlen();
	anchor_mlen = divm[anchor]->getMatchlen();
	nglen = divm[i]->getNogaplen();
	anchor_nglen = divm[anchor]->getNogaplen();
	//tcov = divm[i]->getTcov();
	//anchor_tcov = divm[anchor]->getTcov();

	if ( (zval[i][divstat::idxidentity()] > IZC ||
			zval[i][divstat::idxngidentity()] > IZC ||
			zval[i][divstat::idxsimilarity()] > IZC) &&
			qcov < anchor_qcov) 
	{
		score = divm[i]->getScore() * anchor_qcov/qcov;
		mlen = mlen * anchor_qcov/qcov;
		nglen = nglen * anchor_qcov/qcov;
		if (score > 0.9*divm[anchor]->getScore() ||
				mlen > 0.9*anchor_mlen ||
				nglen > anchor_nglen) 
			return true;  // is due to partial target sequence
		// not 100% right!! it is just a wild guess!!
	}
	return false;
}

bool gdiagnosis::anchorPartial() const {
	pair<int, double> max_qcov = maxqcov();
	pair<int, double> max_iden = maxidentity();
	pair<int, double> max_ngiden = maxngidentity();
	pair<int, double> max_sim = maxsimilarity();

	if ( divm[anchor]->getQcov() < max_qcov.second &&
			divm[anchor]->getIdentity() > 0.9*max_iden.second &&
			divm[anchor]->getNgidentity() > 0.9*max_ngiden.second &&
			divm[anchor]->getSimilarity() > 0.9*max_sim.second) 
	{
		double anchor_qcov = divm[anchor]->getQcov();
		double score, mlen, nglen;
		score = divm[anchor]->getScore() * max_qcov.second/anchor_qcov;
		mlen = divm[anchor]->getMatchlen() * max_qcov.second/anchor_qcov;
		nglen = divm[anchor]->getNogaplen() * max_qcov.second/anchor_qcov;

		// testing corrected features of the anchor division
		if (score > 0.9*maxscore().second || 
				mlen > 0.9*maxmatchlen().second ||
				nglen > 0.9*maxnogaplen().second )
		{
			return true;
		}
	}
	return false;
}

///////////// output function ////////////////////////////////////
ostream& operator<<(ostream& ou, const gdiagnosis &g) {
	operator<<(ou, (group)g);  // it calls the base class operator<<()
	if (g.divCnt < 2) return ou;
	// do not output anything besides the group info
	
	int i,j;
	ou << "stat: (mean,std)\n";
	for (i=0; i<FEAT; i++) ou << divstat::features[i] << '\t';
	ou << endl;
	for (i=0; i<FEAT; i++) 
		ou << setprecision(2) << g.stat[i][0] << "," << g.stat[i][1] 
			<< '\t';
	ou << endl;
	ou << "zval:\n";
	for (i=0; i <= gconst::divisions.size(); i++) {
		if (g.divm[i]) {
			if (i<g.divisions.size()) 
				ou << setw(7) << g.divisions[i] << '\t';
			else ou << setw(7) <<  "sum abs\t";
			for (j=0; j<=FEAT; j++) ou << g.zval[i][j] << '\t';
			ou << endl;
		}
	}
	ou << "histogram:\n";
	for (i=0; i<4; i++) {
		switch (i) {
			case 0:
				ou.width(8);
				ou << "low:\t";
				break;
			case 1:
				ou.width(8);
				ou << "center:\t";
				break;
			case 2:
				ou.width(8);
				ou << "high:\t";
				break;
			case 3:
				ou.width(8);
				ou << "sum:\t";
				break;
			default:
				ou << "grograming error, " << i << " not a possible value\n";
				break;
		}

		for (j=0; j<FEAT; j++) ou << g.hgrm[i][j] << '\t';
		ou << endl;
	}
	return ou;
}

/*  not needed
ostream& gdiagnosis::printRemove(ostream& ou) {
		ou << "removed div: ";
		copy(poordiv.begin(), poordiv.end(), ostream_iterator<int>(ou, " "));
		for (int i=0; i<poordiv.size(); i++) {
			ou << divisions[ poordiv[i] ] << " ";
		}
		ou << endl;
		return ou;
}
*/
// for human inspection and program design. Data on top
// statistics at buttom
ostream& gdiagnosis::dumpWithZval(ostream &ou) const {
	int i;
	for (i=0; i<divisions.size(); i++) 
		if (divm[i]) {
			ou << query << '\t' << *divm[i] << '\t' 
				<< setprecision(3) << setw(6);
			for (int j=0; j<FEAT+1; j++) ou << zval[i][j] << '\t';
			if (i != pivotdividx) {
				if (zval[i][divstat::idxscore()] < -1.2) ou << " <-***\n";
				else if (zval[i][divstat::idxscore()] > 1.2) ou << " ***->\n";
				else ou << endl;
			}
			else ou << endl;
		}
	ou << "summary\tstat:\t";
	for (i=0; i<FEAT; i++) ou << setprecision(1) << stat[i][divstat::idxidentity()] << "|" 
		<< setprecision(2) << stat[i][divstat::idxmatchlen()] << '\t';
	for (i=0; i<FEAT+1; i++) ou << zval[divisions.size()][i] << '\t';
	ou << endl;
	return ou;
}

void gdiagnosis::getspace() {
	zval.resize(divisions.size()+1);
	norm.resize(divisions.size());
	for (int i=0; i< divisions.size(); i++) {
		zval[i].resize(FEAT+1);
		norm[i].resize(FEAT);
	}
	zval[divisions.size()].resize(FEAT+1); // norm and zval diff. dim.
}

//////////////// gstat class ///////////////////////////////
//
///////// initialize static members ///////////////
//int gstat::numdiv=11;  // is it the count of target divisions?
// if the group object has not been initialized, it may cause 
// problem!!
//int gstat::numdiv=group::divisions.size();  // is it the count of target divisions?
//string gstat::pivotDiv="human";
//string gstat::pivotDiv=group::pivotdiv;

/* numdiv number of divisions included in this class 
 * Construct an default object by constructing an 3-d vector
 * sum */
gstat::gstat() : groupCount(0), guid_produced(false) {
	cnt.resize(divisions.size());
	sum.resize(divisions.size());
	guid.resize(divisions.size());
	for (int i=0; i<divisions.size(); i++) {
		cnt[i]=0;
		sum[i].resize(FEAT);
		guid[i].resize(FEAT);
		for (int j=0; j<FEAT; j++) {
			sum[i][j].resize(2, 0.0);  // resize and initialize
			guid[i][j].resize(2, 0.0);
		}
	}
}

/* for training the picker without guid.
 * */
void gstat::accumulate(const group &g) {
	// must have two or more divisions
	if (g.getDivCnt() < 2) return;

	int si = g.getAnchor();  // index for standardization
	double tmp, avg, Mi;
	int i, j;
	groupCount++;
	/* human is normalized to 100.  This will make the calculation 
	 * more accurate and results more understandable. 
	 * */
	if (si == pivotdividx) {
		for (i=0; i<divisions.size(); i++) {
			if (g.divm[i]) {
				cnt[i]++;      // human div is also counted, but not used
				if (i != si) {
					for (j=0; j<FEAT; j++) {
						tmp = g.divm[i]->feat[j]*pivotvalue/g.divm[si]->feat[j];
						// new method
						if (cnt[i] == 1) { // first element
							sum[i][j][0] = tmp;
							sum[i][j][1] = 0;
						}
						else { 
							Mi = sum[i][j][0];
							sum[i][j][0] = Mi + (tmp - Mi)/(cnt[i]+1);
							sum[i][j][1] = cnt[i]*(sum[i][j][1]/(cnt[i]+1) + pow((sum[i][j][0]-Mi),2));
						}
						
						//sum[i][j][0] += tmp;      // computes average
						//sum[i][j][1] += tmp*tmp;  // computes variance
					}
				}
			}
		}
	}
	else {
		for (i=0; i<divisions.size(); i++) {
			// the anchor is not used for statistics!
			if (g.divm[i] && i != si) { 
				cnt[i]++;
				for (j=0; j<FEAT; j++) {
					if (cnt[si] == 0) { 
						avg = pivotvalue;  // kind of like pseudo count
						//exit(1);
					} // debug
					//else avg = sum[si][j][0]/cnt[si];
					else avg = sum[si][j][0];  // with new formula, it is mean already
					tmp = g.divm[i]->feat[j]*avg/g.divm[si]->feat[j];
					if (cnt[i] == 1) { // first element
						sum[i][j][0] = tmp;
						sum[i][j][1] = 0;
					}
					else { 
						Mi = sum[i][j][0];
						sum[i][j][0] = Mi + (tmp - Mi)/(cnt[i]+1);
						sum[i][j][1] = cnt[i]*(sum[i][j][1]/(cnt[i]+1) + pow((sum[i][j][0]-Mi),2));
					}

					//sum[i][j][0] += tmp;
					//sum[i][j][1] += tmp*tmp;
				}
			}
		}
	}
}

/* training with a guid, guid is a member of the class
 * loaded after the first training without the guid
 * */
void gstat::accumulateWithGuid(const group &g) {
	// must have two or more divisions
	if (g.getDivCnt() < 2) return;

	int si = g.getAnchor();  // index for standardization
	double tmp, Mi;
	int i, j;
	groupCount++;
	/* human is normalized to 100.  This will make the calculation 
	 * more accurate and results more understandable. 
	 * */
	if (si == pivotdividx) {
		for (i=0; i<divisions.size(); i++) {
			if (g.divm[i]) {
				cnt[i]++;      // human div is also counted, but not used
				if (i != si) {
					for (j=0; j<FEAT; j++) {
						tmp = g.divm[i]->feat[j]*pivotvalue/g.divm[si]->feat[j];
						if (cnt[i] == 1) { // first element
							sum[i][j][0] = tmp;
							sum[i][j][1] = 0;
						}
						else { 
							Mi = sum[i][j][0];
							sum[i][j][0] = Mi + (tmp - Mi)/(cnt[i]+1);
							sum[i][j][1] = cnt[i]*(sum[i][j][1]/(cnt[i]+1) + pow((sum[i][j][0]-Mi),2));
						}
						//sum[i][j][0] += tmp;      // computes average
						//sum[i][j][1] += tmp*tmp;  // computes variance
					}
				}
			}
		}
	}
	else {
		for (i=0; i<divisions.size(); i++) {
			// the anchor is not used for statistics!
			if (g.divm[i] && i != si) { 
				cnt[i]++;
				for (j=0; j<FEAT; j++) {
					tmp = g.divm[i]->feat[j]*guid[i][j][0]/g.divm[si]->feat[j];
					if (cnt[i] == 1) { // first element
						sum[i][j][0] = tmp;
						sum[i][j][1] = 0;
					}
					else { 
						Mi = sum[i][j][0];
						sum[i][j][0] = Mi + (tmp - Mi)/(cnt[i]+1);
						sum[i][j][1] = cnt[i]*(sum[i][j][1]/(cnt[i]+1) + pow((sum[i][j][0]-Mi),2));
					}
					//sum[i][j][0] += tmp;
					//sum[i][j][1] += tmp*tmp;
				}
			}
		}
	}
}

/* make sure the calguid() is called before calling this */
ostream& operator<<(ostream &o, gstat &gs) {
	if (!gs.guid_produced) gs.calguid();  // incase it has not been calculated
	o.setf(ios::showpoint);
	o.precision(3);
	o.setf(ios::fixed, ios::floatfield);
	o << "#" << gs.groupCount << " query groups used\n#feature\t";
	int i, j;
	for (i=0; i<FEAT; i++) {
		o << divstat::features[i] << "\t\t";
	}
	o << "\n#div\t";
	for (i=0; i<FEAT; i++) {
		o << "avg\tstd" << "\t";
	}
	o << "n\n";

	//o << "avg\tstd\t for five features  n\n";
	double avg, var, std;
	for (i=0; i<gs.divisions.size(); i++) {
		o << group::divisions[i] << '\t';
		for (j=0; j<FEAT; j++) {
			o << gs.guid[i][j][0] << '\t' << gs.guid[i][j][1] << '\t';
		}
		o << gs.cnt[i] << endl;
	}
	return o;
}
/* calculates the guid matrix from the sum matrix */
double gstat::calguid() {
	int i,j;
	// computes the std from variance
	for (i=0; i<divisions.size(); i++) {
		for (j=0; j<FEAT; j++) {
			if (sum[i][j][1] > 0) sum[i][j][1] = sqrt(sum[i][j][1]);
		}
	}
	// now the sum matrix contains the mean,std of each feature from all divisions
	
	double diff=0;  // difference of mean + std of next round training from first round

	for (i=0; i<divisions.size(); i++) {
		if (i==pivotdividx) {  // no need for computation
			for (j=0; j<FEAT; j++) {
				guid[i][j][0] = pivotvalue;
				guid[i][j][1] = 0;
			}
		}
		else {
			for (j=0; j<FEAT; j++) {
				if (cnt[i] == 0) {  // no value available for this div
					guid[i][j][0] = pivotvalue;
					//cerr << csvStdHighVar << " the csvStdHighVar value\n";
					guid[i][j][1] = csvStdHighVar*2;  // estimate
				}
				else {
					//guid[i][j][0] = sum[i][j][0]/cnt[i];
					//guid[i][j][1] = sqrt(sum[i][j][1]/cnt[i] - pow(guid[i][j][0], 2));
					diff += fabs((sum[i][j][0] - guid[i][j][0]));  // mean diff
					guid[i][j][0] = sum[i][j][0];
					diff += fabs((sum[i][j][1] - guid[i][j][1]));  // std diff
					guid[i][j][1] = sum[i][j][1];
				}
			}
		}
	}
	guid_produced = true;
	return diff;
}

void gstat::zeroguid() {
	for (int i=0; i<divisions.size(); i++)
		for (int j=0; j<FEAT; j++) 
			for (int k=0; k<2; k++) guid[i][j][k]=0;
}
void gstat::zerosum() {
	for (int i=0; i<divisions.size(); i++) {
		cnt[i]=0;
		for (int j=0; j<FEAT; j++) 
			for (int k=0; k<2; k++) sum[i][j][k]=0;
	}
	guid_produced = false;
}
//////////////////////////////////////////////////////////////
///////////////// the trainer class //////////////////////////
//////////////////////////////////////////////////////////////
//trainer::trainer(int numdiv, const string &pivotdiv) {
trainer::trainer() {
	//gstat::setNumdiv(numdiv);  // assignment too late
	// gstat object already constructed
	//gstat::setPivotDiv(pivotdiv);
}

/* the training should be a class method
 * train the picker with training data from well annotated sets */
//bool trainer::train(const string &inf, const string &guidfile) {
bool trainer::train(const string &inf) {
	ifstream IN(inf.c_str());
	if (IN.fail()) { 
		cerr << "cannot open " << inf << "\n";
		return false;  // train failed
	}
	cerr << " reading from file: " << inf << endl;
	ofstream OU("better.tab");  // purify result
	//ofstream ED("marked.tab");  // how it was purified
	//gstat model, conservedModel;
	//
	/* We keep two lists of groups 
	 * for repeated training */
	vector<group> conserved, normal;
	
	//cerr << "Reading conf file ...\n";
	//group::readconf(guidfile);
	//gconst::readconf(guidfile);

	int query;
	//double SMratio;
	IN >> query;
	try {
		while (1) {
			group ag(query, IN);
			if (ag.getDivCnt() < 2) continue;
			//cout << ag << endl;
			//ED << ag << endl;
			//ag.discardLow(ED);
			ag.dumpAsTable(OU);
			//SMratio =  ag.getScoreSMratio(); 
			OU << setprecision(4); 
			ag.dumpAllSMratio(OU);
			if (ag.isConserved()) {
				conserved.push_back(ag);
				OU << " **Conserved ";
				csvdmodel.accumulate(ag);
			}
			else {
				model.accumulate(ag);
				normal.push_back(ag);
			}
			OU << endl;
			OU << endl;
		}
	}
	catch (inputend) {
		cerr << "input ended\n";
		// is the last group added to the statistics
		//return 0;
	}
	IN.close();
	OU.close();

	/////////  output model //////////////////////
	model.calguid();
	csvdmodel.calguid();
	//cout << "# Model from training without guid, first round\n"
	//		<< "############ Regular Model ###############\n"
	//		<< model << endl;
	//cout << "######  Conserved Model #############\n"
	//		<< csvdmodel << endl;
	OU.open("picker.guid");
	OU << "# Training from data without guid\n# Regular Model\n" << model << endl;
	OU << "# Conserved Model\n" << csvdmodel << endl;

	cerr << " beter.tab contains purified results\n"
		<< " marked.tab contains marked results\n";

	//// continued training //////////////////////
	int i, k;
	double diff1=1000, diff2=1000;  // give it a large number to start the loop
	for (k=0; (diff1>0.02 or diff2>0.2) && k<100; k++) {  // repeat twice
		model.zerosum();
		for (i=0; i<normal.size(); i++) {
			model.accumulateWithGuid(normal[i]);
		}
		diff1 =  model.calguid();
		//cout << "# Diff=" << diff << "####### model after " 
		//	<< k+2 << " rounds\n" << model << endl;
		//OU << "####### Model after " << k+2 << " rounds, Diff=" << diff1 << "  "
		//	<< model << endl;

		csvdmodel.zerosum();
		for (i=0; i<conserved.size(); i++) {
			csvdmodel.accumulateWithGuid(conserved[i]);
		}
		diff2 = csvdmodel.calguid();
		//cout << "# Diff=" << diff << "### conserved model\n";
		//cout << csvdmodel << endl << "## ===================================\n";
		//OU << "# Diff=" << diff2 << "  ## Conserved model\n";
		//OU << csvdmodel << endl << "## ===================================\n\n";
	}
	OU << "# Model with diff=" << diff1 << " after " << k+2 << " rounds of training\n"
		<< model << endl;
	OU << "# Conserved Model  diff=" << diff2 << "\n";
	OU << csvdmodel << endl << "## ===================================\n\n";
	cerr << "Model guid written to file picker.guid\n";
	OU.close();
	return true;
}
