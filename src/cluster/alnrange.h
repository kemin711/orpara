#ifndef ALNRANGE_H
#define ALNRANGE_H

#include <iostream>
#include <string>
#include <vector>
#include <list>

#ifdef HAVE_PG
#include "libpq++.h"
#endif

using namespace std;

namespace orpara {

// detect chimera sequence given a cluster number

/** alnrange class
 * the range is an individual, not a collection or average
 * need another object to keep this information
 *
 * Using another data type such as PgDatabase may not be good.
 *
 * Input and output should use standard lib containners
 * instead of database connection object.  This will make
 * the class more portable.
 * 
 * There is a much better and more generic implementation
 * of ragne in range directory. I should rename the range
 * in this directory.
 *
 * Should be alnrange
 *
 * */

class alnrange {
	public:
		alnrange() : begin(0), end(0), score(0), ngidentity(0), cov(0) {}
		alnrange(const alnrange &r);
        /**
         * b: begin of the match; e: end of the match
         * s: score; ng: ngidentity, c: coverage
         */

		alnrange(int b, int e, double s, double ng, double c) 
			: begin(b), end(e), score(s), ngidentity(ng), cov(c) {}

#ifdef HAVE_PG
        /** special for PgDatabase */
		alnrange(PgDatabase* db, int i);
#endif

		virtual ~alnrange() {}
		alnrange& operator=(const alnrange& r);

		/* this range overlap with another one
		// This methods have defects, and should be applied
		// on large range (>margin) only
		*/
		bool overlap(const alnrange &r, const int margin);

		/* use internal class wide parameters */
		bool overlap(const alnrange &r);
		int length() const { return end-begin+1; }
		int getBegin() const { return begin; }
		int getEnd() const { return end; }
		double getScore() const { return score; }
		double getNg() const { return ngidentity; }
		double getCov() const { return cov; }

		/* out put all fields
		 * begin,end,ngidentity,cov; separator: "\t"
		 * */
		virtual string asTabedString() const;
		virtual string asDelimitedString(const char sep[]=",") const;

		static string fields();

		/* pro produce only the essential information for debug
		 * tab delimited string output
		 * begin,end,ngidentity
		 * */
		virtual string essentialInfo() const;

		friend ostream& operator<<(ostream &ous, const alnrange& r);

		static void setmargin(int olpcut, float olpfrac) {
			ovlpcut=olpcut; ovlpfraction=olpfrac; }

	protected:
		int begin; // match begin index, usually 1-based
		int end;  // matche end index
		// keeps average score, 
      /** addable score from blast or alignment programs */
		double score;
      /** non-gapped identity, identity should also be fine */
		double ngidentity;
      /** coverage on query */
		double cov; // coverage relative to query

        static int ovlpcut;
        static float ovlpfraction;
	//double identity;
};

/** did not implement the copy constructor
 * we will use the default bitwise supplied by the
 * compiler
 *
 * This simply a match two ranges. Could be implemented better.
 */
class rangePair : public alnrange {
	public:
		rangePair() : alnrange(), tprtid(), tcov(0), tmodelid(0),
			tgenomicid(), tstrand(), tstart(-1), tend(-1) {}
		rangePair(const rangePair& r);
		/**
		 * tpi: target proteinid
		 * tc: target coverage
		 * tmi: target model id
		 * tgi: target genomic id
		 * tst: target strand
		 * ts: target start, te: target end
		 */
		rangePair (int b, int e, double s, double ng, 
				double c, const string &tpi, double tc, 
				int tmi, const string &tgi, char tst, 
				int ts, int te) 
			: alnrange(b,e,s,ng,c), tprtid(tpi), tcov(tc),
			tmodelid(tmi), tgenomicid(tgi), tstrand(tst), 
			tstart(ts), tend(te) {}
		double getTcov() const { return tcov; }
		string getTprtid() const { return tprtid; }
		/**
		 * Out all fields in tab delimited format
		 * begin,end,ngidentity,cov,
		 * tprtid, tmodelid, tgenomicid, tstrand,
		 * tstart, tend, tcov
		 * */
		string asTabedString() const;
		string asDelimitedString(const char sep[]=",") const;
		string genomicInfo(const char sep[]="\t") const { 
			return tgenomicid + sep + tstrand; }
		string getGenomic() const { return tgenomicid; }
		char getStrand() const { return tstrand; }
		static string fields();

		/** output string with field separator TAB '\t'
		 * begin, end, score, ngidentity, cov, tprtid, tmodelid, tstart, tend
		 * */
		string essentialInfo() const;

		~rangePair() {}

		bool sameGene(const rangePair &r) const;
		/**
		 * return a list of all fields.
		 * useful for generating database tables
		 * not very useful, itoa not implemented, too much trouble
		 * to implement this function 
		 * */
		//list<string> fieldsAsList() const;

		// before using this class, this parameter must be
		// set to a proper value
		static void setIntronLimit(int length) { distance_cut=length; }

	private:
		/* add target spcific staff */
		string tprtid; // target id
		double tcov; // target coverage
		int tmodelid;
		string tgenomicid;
		char tstrand; // + or -
		int tstart;  // genomic start
		int tend;   // genomic end

		static float ngdiff_cut;

		/** This parameter is related to the mean / stddev of the 
		 * intron size from this particular organism
		 **/
		static int distance_cut;
};

class SplitResult {
	private:
		string guide;
		int guideLen;
		list< pair<rangePair,rangePair> > joins;

	public:
		SplitResult() {}
		void setGuide(const string &g, int len) { guide=g; guideLen=len; }
		void add(const rangePair *left, const rangePair *right) {
			joins.push_back(make_pair(*left, *right)); }
		bool empty() const { return joins.empty(); }
		// for constructing SQL use sep=",", for table dump use \t
		list<string> outputRow(const char sep[]=",") const;
};

// begin and end become the outer-most value
/** class avgrange
 * This class have design problem.  I have put the parameters
 * into the range.  There should be a hierchy of classes
 * to derive from the base class
 */
class avgrange {
	public:
		// default constructor, build an empty object
		avgrange() 
			: begin(999999), end(0), 
           sumbegin(0), sumend(0),
           sumscore(0), sumng(0), sumcov(0), 
           n(0), members(), 
			  covs(), sorted(false) { }

		/* construct an avgrange object out of r
		 * avgrang is the same as range if there is only one member.
		 * It makes a copy from r.
		 * r must be created by the new operator.
		 * */
		avgrange(const alnrange &r);
		avgrange(const alnrange *rp);
		~avgrange();
		bool overlap(const alnrange &r, const int margin=10);

		// this is the most useful method, for accumulating
		// overlapping ranges
		//void merge(const range &r);
		void merge(const alnrange *r);
		/* do pointer manipulation, no object copying
		 */
		void merge(const avgrange *ar);

		double length() const { return (sumend-sumbegin+n)/static_cast<double>(n); }
		double getBegin() const { return sumbegin/static_cast<double>(n); }
		double getEnd() const { return sumend/static_cast<double>(n); }
		double getScore() const { return sumscore/sumcov; }
		double getNg() const { return sumng/sumcov; }
		double getCov() const { return sumcov/n; }
		int maxlength() const { return end-begin+1; }
		int minbegin() const { return begin; }
		int maxend() const { return end; }
		int getCount() const { return members.size(); }
		const vector<const alnrange* > & getMembers() const { return members; }
		// for human to read
		friend ostream& operator<<(ostream &ous, const avgrange &ar);

        /* output one line in table format, not line terminator, 
		 * The fields are defined in the colheaders() method
		 **/
		ostream& writeTable(ostream &ous) const;
		/** This method is an analogue of the above method.
		 * It returns a tab-delimited string
		 * */
		string asDelimitedString(const char sep[]=",") const;

		/* with SQL comment format, out put range information
		 * */
		ostream& sqlinfo(ostream &ous);

		// assume covs is sorted
		double getMedianCov() const;

		/* at testing stage write output to stdout
		 * This function test all pair-wise split genes, 
		 * 3 or more split is rare, but is present can be 
		 * easily detected. for example
		 * A-B B-C will be combined to A-B-C
		 * When we say it is a split gene, you have to make
		 * sure that the standard is not a chimera!
		 * */
		string checkSplit_debug(const avgrange &r) const;
		list<string> checkSplit(const avgrange &r, const char sep[]=",") const;
		SplitResult testSplit(const avgrange &r, const char sep[]=",") const;

		// for chimera detection
        static string colheaders() {
            return "min_begin\tmax_end\tavg_begin\tavg_end\tn\tavg_score\tavg_ngidentity\tavg_coverage\tmedian_coverage";
        }

	private:
		int begin, end; // outer most value
		int sumbegin, sumend;
		double sumscore, sumng, sumcov;
		//double sumweight;
		int n;
		/** all the ranges in this collection.  No particular
		 * ordering.
		 * */
		vector<const alnrange*> members;

      // caching variables
		mutable vector<double> covs; 
		/* used to compute median coverage
		 * Only after sorting, median can be computed
		 * This is a state tag for covs
		 */
		mutable bool sorted;

		//int seqlength;  // length fo the complete sequence
		//string id;  // sequence identifier
		// as oppose to matchlen
};
}
#endif
