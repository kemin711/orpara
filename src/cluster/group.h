#ifndef GROUP_H
#define GROUP_H
/** version 2.0 for discriminating against divisions that don't have
 * a true match */
//#define HAVE_NAMESPACE_STD

#include <iostream>
#include <cstring>
#include <utility>
//#include <libpq++.h>  // for database input
#include <pqxx>  // for database input
// change to the libpqxx interface
#include <string>
#include <vector>
#include <fstream>
#include "gconst.h"

using namespace pqxx;

//#define DIV 11    // fr to 11 other dbs comparison; this is stored
//inside divisions vector
#define FEAT 8  // we use 8 features to characterize each division

class inputend {};  // for signaling

/** A simple statistical representation of the max(...) features from
 * each division.  In the first test run I used 5 features.  Now I have
 * a much richer features:
 * max of 
 * 	 0  qcov float4,
 * 	 1  tcov float4,
 * 	 2  score float,
 * 	 3  identity float4,
 * 	 4  ngidentity float4,
 * 	 5  matchlen integer,
 * 	 6  nogaplen integer,
 * 	 7  similarity float4,
 * 	    ======
 * all the features except the simple feature that should be count
 * count(distinct(simple)) if == 2 means mixed and lower quality of
 * some of the hits.
 *
 * This class simply holds the values.
 * */
class divstat
{
	public:
		divstat();  // default constructor

		/* Construct object from from a stream. Input pointer read to read 
		 * div, qcov, tcov, score, identity, ngidentity, similarity, 
		 * matchlen, and nogaplen.
		 * All caller should use the interface functions to use
		 * this object.
		 * */
		divstat(istream &in);  // use a get function
		void read(istream &in);

		/** constructor from database input
       * now the result type of the pqxx
       * */
		//divstat(PgDatabase &pgdb, int row);
		divstat(result &qres, int row);

		/* copy constructor */
		divstat(const divstat &d);

		string getid() { return string(id); }

		/* testing the div name */
		bool operator==(const char s[]) { return !strcmp(s, id); }
		friend ostream& operator<<(ostream &o, const divstat &ds);
		divstat& operator=(const divstat &d);
		bool operator>(const divstat &ds) const;
		bool operator<(const divstat &ds) const;

		/* this object is smaller than d after d is multiplied by k.
		 * if d is 0 or null, then returns false. */
		bool smallerThan(const divstat &d, double k) const;
		bool smallerThan(const divstat *const d, const double k) const;
		void scale(double f);
		/* use the following functions for flexibility */
		double getQcov() const { return feat[0]; }
		double getTcov() const { return feat[1]; }
		double getScore() const { return feat[2]; }
		double getIdentity() const { return feat[3]; }
		double getNgidentity() const { return feat[4]; }
		double getSimilarity() const { return feat[5]; }
		double getMatchlen() const { return feat[6]; }
		double getNogaplen() const { return feat[7]; }
		ostream& dump(ostream &ou);
		/* the following functions return the index of features in the 
		 * features array
		 * */
		static int idxqcov() { return 0; }
		static int idxtcov() { return 1; }
		static int idxscore() { return 2; }
		static int idxidentity() { return 3; }
		static int idxngidentity() { return 4; }
		static int idxsimilarity() { return 5; }
		static int idxmatchlen() { return 6; }
		static int idxnogaplen() { return 7; }

	//private:
		char id[5];         // the target div name such as hs, bony, etc.
		double feat[FEAT];  // feature values, fixed according to the above order
		static const char* features[FEAT];  // feature string
};

/* A class to represent one query match against all target divisions
 * This program needs the guid file to instruct what is the query div
 * and what are target divisions. 
 *
 * The first action of this class is to read and initialize data from the
 * guid file*/
class group : public gconst
{
	public:
		friend class gstat;

		/* default constructor, allocate space, initialize to impossible
		 * values*/
		group();
		group(const group &g);
		group& operator=(const group &g);
		~group() { for(int i=0; i<divisions.size(); i++) delete divm[i]; }

		//this way of noticing the end of input stopped working 
		group(int &q, istream &in) throw(inputend);
		// replacement for special constructor, return true if more
		// objects in the input stream, q is the query hint
		bool next(int &q, istream &in);

		////////////  Output Functions  //////////////////////
		
		/* for human to read */
		friend ostream& operator<<(ostream &o, const group &g);

		/* output the raw data as tab delimited, for database entry */
		void dumpAsTable(ostream &ou);

		/**
		 * Only output the (query,div) pair, no feature values.
		 * For database loading.
		 */
		void dumpKey(ostream &ou);

		// dump everything to ou for computer to read
		ostream& dump(ostream &ou);
		// output the stat matrix
		void dumpStat(ostream &ou);
		void dumpAllSMratio(ostream &ou);

		/////////  Access Functions ///////////////
		// the first is index, second is the value of max tcov
		pair<int, double> maxqcov() const { return max(divstat::idxqcov()); }
		pair<int, double> maxtcov() const { return max(divstat::idxtcov()); }
		pair<int, double> mintcov() const;
		pair<int, double> maxidentity() const { return max(divstat::idxidentity()); }
		pair<int, double> maxngidentity() const { return max(divstat::idxngidentity()); }
		pair<int, double> maxsimilarity() const { return max(divstat::idxsimilarity()); }
		pair<int, double> maxmatchlen() const { return max(divstat::idxmatchlen()); }
		pair<int, double> maxnogaplen() const { return max(divstat::idxnogaplen()); }
		pair<int, double> maxscore() const { return max(divstat::idxscore()); }
		pair<int, double> max(int f) const;

		/* obtain the minimum target sequence coverage, set the division index
		 * in the divm array.  */
		double getMinTcov(int &div) const;

		int getDivCnt() const { return divCnt; }
		// obtain the current anchor index
		int getAnchor() const { return anchor; }
		// obtain the next anchor, returns -1 if no more anchor available
		// for example, the first anchor human is bad, we want to get the 
		// next one. If we have only two divisions, and we want the next again
		// we will be running out of anchors!
		// s is the starting poing for anchor
		int nextAnchor();

		/* SM means standard deviation/mean ratio */
		double getIdenSMratio() const { return stat[divstat::idxidentity()][1]/stat[divstat::idxidentity()][0]; }
		double getMlenSMratio() const { return stat[divstat::idxmatchlen()][1]/stat[divstat::idxmatchlen()][0]; }
		double getQcovSMratio() const { return stat[divstat::idxqcov()][1]/stat[divstat::idxqcov()][0]; }
		double getTcovSMratio() const { return stat[divstat::idxtcov()][1]/stat[divstat::idxtcov()][0]; }
		double getScoreSMratio() const { return stat[divstat::idxscore()][1]/stat[divstat::idxscore()][0]; }

		/* If identities are higher than 0.85 and SMratio < 0.1 or SMratio < 0.2 for
		 * bony and vrt 
		 * as defined in *.istr file under
		 * CONSERVED 0.85 0.1 0.2 bony vrt
		 * */
		bool isConserved() const;

		////// Mutating Functions ///////////////////////////////////

		/* Remove a target division identified by index position i
		 * from the group output information to ou for debuging. 
		 * Recalculate stat. If anchor removed, reposition anchor.
		 * */
		void rmbranch(int i, ostream &ou); // debug version
		void rmbranch(int i);
		
	protected:
		bool highVarDivPresent() const;
		/* this function compares divisions i with division b to e
		 * assuming i is mostly related.   This one is not very useful
		 * because it is based on human as query and the targets are
		 * in order mostly related to human */
		bool checkAndRemove(int i, int b, int e, ostream &ou);

		/* initialize or recalculate stat[5][2]; call this one after
		 * the object has been changed by removing a branch of division.
		 *
		 * Stattistics of each feature over all existing matching target
		 * divisions
		 * */
		void dostat();

		int query;  // the cds_id of a particular query
		/* the divm array has a fixed order as defined in the guid file
		 * It hold the max values from each matching division
		 * */
		vector<divstat*> divm;
		int divCnt;  // the number of target divisions (including anchor) for this query
		int anchor;       //index for division picked for normalization 
		/* (avg,std) of features over divisions in each group */
		double stat[FEAT][2]; 
};

/** for group diagnosis of quality control 
 * this class will work on individual groups of divisions and
 * figure out whether they are the same orthologue cluster or not
 * */
class gdiagnosis : public group {
	public:
		/* Needs to allocate space for copying, zcut set to 1 */
		gdiagnosis();
		gdiagnosis(double cut);
		gdiagnosis(const gdiagnosis &gd);
		// construct objects from in
		//gdiagnosis(int &q, istream &in) throw(inputend);
		
		// for examining out put in the context of statistics
		ostream&  dumpWithZval(ostream &ou) const;
		//void normalize();  // use internally stored guid
		bool next(int &q, istream &in);

		void setzcut(double c) { zcut=c; }
		double getzcut() const { return zcut; }

		///////// output function //////////////////////

		gdiagnosis& operator=(const gdiagnosis& g);
		friend ostream& operator<<(ostream &ou, const gdiagnosis &g);
		//ostream& printRemove(ostream &ou);
		
		/////////////  Diagnostic Functions ///////////////////
		// the most vigorous test, all features z-values are smaller
		// than zcut

		// there is no low and high != sum (anchor should not be low)
		bool qcovpass() { return hgrm[0][divstat::idxqcov()] == 0 && hgrm[2][divstat::idxqcov()] != hgrm[3][divstat::idxqcov()]; }
		bool scorepass() { return hgrm[0][divstat::idxscore()] == 0 && hgrm[2][divstat::idxscore()] != hgrm[3][divstat::idxscore()]; }
		// any of the identity passed, and the anchor should also pass
		bool qualitypass();

		bool passed(const double zcut) const;
		/* divCnt == 1 is good; decision should be made outside this class
		 * test divCnt == 1 or conserved or passedAvg
		 * */
		bool goodQuality(double zcut) const;

		bool passedAvg(const double zcut) const;
		bool passedQcov(const double zcut) const; // only test qcov
		/* check all three identity (identity, nogap, similarity), if anyone
		 * passed then passed */
		bool passedIdentity(const double zcut) const;
		/* check for nogaplength, or matchlength, if anyone passed then passed */
		bool passedLength(const double zcut) const;

		bool passedall(const double zcut) const {
			return (passed(zcut) || goodQuality(zcut) || passedAvg(zcut) ||
				passedQcov(zcut) || passedIdentity(zcut) || passedLength(zcut));
		}

		/* some of the low scores is due to partial sequence, that we name
		 * coverageDefect */
		bool coverageDefect(double zcut) const;
		bool lastTest() const;

		/////////////// helper functions /////////////////////////
		void calzval();   // helper function, initializes zval and norm

		// remove divisions that are low, recalculate all derived matrices
		// return a list of removed divisions (index value in divm vector)
		vector<int> rmlow(double lowercut=-3, double avgzc=3); // check the quality of each division, load remove array

		/* return the average of all zvalues. It is signed. 
		 **/
		double getAvgZval() const { return zval[divisions.size()][FEAT]/(FEAT*(divCnt-1));}
		/* same as above just easier to type */
		double avgzval() const { return zval[divisions.size()][FEAT]/(FEAT*(divCnt-1));}

		bool targetPartial(int i) const;  // division i has partial sequence
		bool anchorPartial() const;
		
		//bool improve(ostream &ou);
		bool trimAndTest(ostream &ou, double zcut);

		void trimByIden(); // this method cannot be applied to 
		// other situations than hs->fr hs->non-fish divs

		/* d if for debug.  It contains the fixed object */
		bool fixScoreAndTest(double szcut, gdiagnosis &d) const;

		/////////////// static section //////////////////////////////
		static void readGuid(istream &in);   // needs to be removed
		//static double guid[8][5][2];
		/* the guid is the result from trainer class */

	private:
		static vector< vector< vector<double> > > guid;
		/* standard normal z value as compared to the guid statistics
		 * dimension: (numdiv+1, numfeat+1) or (divisions.size()+1, FEAT+1) 
		 * Calculated after construction from input stream.
		 * the bottom row is the sum(fabs(z))/divCnt-1
		 * the last column is the sum of fabs(z) of all features for each div
		 * Some rows may not be populated (missing divisions).
		 * */
		vector< vector<double> > zval; //division.size()+1 rows, FEAT+1 columns  

		/* normalized features of populated divisions, with divisions.size()
		 * row, and FEAT columns */
		vector< vector<double> > norm;

		int hgrm[4][FEAT];  // histogram of low, passed, high, and sum counts, 
		double zcut;        // blow this value is considered good
		//vector<bool> remove;    // division is removed or not
		//vector<int> poordiv;  // list of div needs to be removed
		void getspace(); // helper function to allocate space for zval
	                    //	and norm, all elements has value 0
};

/** produces the guid matrix from training input data 
 * the total number of divisions should be included in constructing
 * each object
 * */
class gstat : public gconst {
	public:
		gstat();

		/* for training the picker only
		 * training is done with some averaging of previous data points
		 * this may cause big errors.  Repeated training may be implemented
		 * in the future to increase the accuracy of the training. 
		 * */
		void accumulate(const group &g);
		/* for repeated training from the same data 
		 * REMEMBER to clear the matrices and cnt vector before calling this
		 * function after either calling accumulate or accumulateWithGuid
		 * */
		void accumulateWithGuid(const group &g);

		/* copy the guid matrix from the sum matrix that already contain
		 * the mean,std pair
		 * will overwrite the old values in guid  matrix
		 * return the sum of differences of old guid and new guid 
		 * matrix
		 * */
		double calguid();
		/* zero the sum matrix and the cnt vector for next round of training
		 * */
		void zerosum();
		/* It seem not necessary to call zeroguid */
		void zeroguid(); // call this for a new round of training
		/* dump the simple statics to a file for use as guid 
		 * if gs have not produced the guid matrix this function will
		 * do the calculation and thus modifying the gs object.*/
		// not needed zerosum does the same void clear();
		friend ostream& operator<<(ostream &o, gstat &gs);
		bool guid_produced;  // the guid matrix is produced or not

		//static void setNumdiv(int nd) { numdiv=nd; }
		static int getNumdiv() { return divisions.size(); }
		//static void setPivotDiv(string pdv) { pivotDiv=pdv; }

		static const int pivotvalue = 10;  // smaller the better

	private: 
		/* add info from normalizing g by Di */
		//bool suck(const group &g, int Di);
		vector<int> cnt;         // for each bin of divisions
		//double sum[8][5][2];   // for calculating mean and variance
		//double*** sum;   // 3-d array for calculating mean and variance
		/* A lazy version of 3D array for 
		 * division, feature, (mean, standard deviation)
		 * this is used for computing the mean and std
		 * */
		vector< vector< vector<double> > > sum;
		// my coding time
		//double sumsqr[8][5];
		vector< vector< vector<double> > > guid;
		int groupCount;         // total groups used

		///////////////    static members  /////////////////////////////////
		//static int numdiv;      // total number of target divisions in input data
//		static string pivotDiv; // the div whose value is set to 100
		// this value is redundant with the group::pivotdiv 
		// need to remove this member in the future
};

/* pure training, will not discard low uising the crude method 
 * If the training dataset is not large, you can add a complete
 * set of training data to the beginning of the training data.  This
 * way, every division will be represented.  For example, there are
 * very little sequence information available for the cartilaginous
 * fish.  */
class trainer {
	public:
		/* numdiv is the number of target divisions
		 * pivotdiv is the pivot division used for normalization
		 * they have the exace meaning as defined in gstat class.  
		 * They are used to set the static member of the gstat class
		 * */
		//trainer(int numdiv, const string &pivotdiv);
		trainer();
		/* inf contains the training data, 
		 * where * is the query division 
		 * the training result will be written to the STDOUT and file
		 * picker.guid
		 * */
		//bool train(const string &inf, const string &guidfile);
		bool train(const string &inf);
	private:
		gstat model, csvdmodel;
};

#endif
