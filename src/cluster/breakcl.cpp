/* File: breakcl.cpp
 *
 * This program use the cluster and self_aln information from two tables
 * stored in the postgres database.  Then it tries to find the weakest 
 * link and remove it and repeat the clustering process.
 *
 * Here I defined Relrow that contain more fields than the ngidentity
 * that is used for clustering.  These extra fields will allow future
 * modifications of this program to work with more complicated methods
 * of clustering.
 * */

//#include "forest.h"
/*
#include <unistd.h>
#include <dirent.h>
#include "scorepair.h"
#include <iterator>
*/

#include <iostream>
#include <cstring>
//#include <libpq++.h>
#include <pqxx>
#include <string>
#include <list>
#include <map>
#include <functional>
#include <algorithm>
#include <set>
#include "stddev.h"
#include "hatrees.h"
#include "strformat.h"
#include <sstream>

using namespace std;
using namespace pqxx;

// the query must be structured as the following query
class Relrow {
	public:
		string query;
		string target;
	 	float qcov;
		float tcov;
		int matchlen;
		float identity;
		float ngidentity; 
		float similarity;
		float score;
	
		Relrow(result &qres, int row);
		Relrow() : query(), target(), qcov(0), tcov(0), matchlen(0), identity(0), ngidentity(0), similarity(0), score(0) {}
		float smcov() const { return qcov<tcov? qcov : tcov; }
		float lgcov() const { return qcov>tcov? qcov : tcov; }
};

Relrow::Relrow(result &qres, int row) {
//	int fields = db.Fields();
//	must have 9 fields
//	   0      1      2       3     4          5         6            7        8
	//query, qcov, target, tcov, matchlen, identity, ngidentity, similarity, score 
   /*
	query = db.GetValue(row, 0);
	target = db.GetValue(row, 2);
	qcov = atof(db.GetValue(row, 1));
	tcov = atof(db.GetValue(row, 3));
	matchlen = atoi(db.GetValue(row, 4));
	identity = atof(db.GetValue(row, 5));
	ngidentity = atof(db.GetValue(row, 6));
	similarity = atof(db.GetValue(row, 7));
	score = atof(db.GetValue(row, 8));
   */

   qres[row][0].to(query);
   qres[row][2].to(target);
   qres[row][1].to(qcov);
   qres[row][3].to(tcov);
   qres[row][4].to(matchlen);
   qres[row][5].to(identity);
   qres[row][6].to(ngidentity);
   qres[row][7].to(similarity);
   qres[row][8].to(score);
}
///////////////////////// end of Relrow class ////////
// predicates
class Relrow_less : public binary_function<Relrow, Relrow, bool> {
	public:
	bool operator()(const Relrow& r1, const Relrow& r2) const {
		return r1.ngidentity < r2.ngidentity; }
};
class Relrow_lessScore: public binary_function<Relrow, Relrow, bool> {
	public: 
		bool operator()(const Relrow& r1, const Relrow& r2) const {
			return r1.score < r2.score; }
};

// smaller coverage
class Relrow_lessSmcov: public binary_function<Relrow, Relrow, bool> {
	public: 
		bool operator()(const Relrow& r1, const Relrow& r2) const {
			return r1.smcov() < r2.smcov(); }
};

// larger coverage
class Relrow_lessLgcov: public binary_function<Relrow, Relrow, bool> {
	public: 
		bool operator()(const Relrow& r1, const Relrow& r2) const {
			return r1.lgcov() < r2.lgcov(); }
};
class Relrow_lessMatchlen: public binary_function<Relrow, Relrow, bool> {
	public: 
		bool operator()(const Relrow& r1, const Relrow& r2) const {
			return r1.matchlen < r2.matchlen; }
};

class Relrow_eq2 : public binary_function<Relrow, Relrow, bool> {
	public:
	bool operator()(const Relrow& r1, const Relrow& r2) const {
		return r1.ngidentity == r2.ngidentity; 
	}
};

class Relrow_eq1 : public unary_function<Relrow, bool> {
	float cutoff;
	public:
	explicit Relrow_eq1(const float cut) : cutoff(cut) { }
	bool operator()(const Relrow& r) const {
		return r.ngidentity == cutoff; 
	}
};

class Relrow_eq1Score : public unary_function<Relrow, bool> {
	private:
		float cutoff;
	public:
		explicit Relrow_eq1Score(const float cut) : cutoff(cut) { }
		bool operator()(const Relrow& r) const {
			return r.score == cutoff; 
		}
};

class Relrow_eq1Smcov : public unary_function<Relrow, bool> {
	private:
		float cutoff;
	public:
		explicit Relrow_eq1Smcov(const float cut) : cutoff(cut) { }
		bool operator()(const Relrow& r) const {
			return r.smcov() == cutoff; 
		}
};
class Relrow_eq1Lgcov : public unary_function<Relrow, bool> {
	private:
		float cutoff;
	public:
		explicit Relrow_eq1Lgcov(const float cut) : cutoff(cut) { }
		bool operator()(const Relrow& r) const {
			return r.lgcov() == cutoff; 
		}
};
class Relrow_eq1Matchlen : public unary_function<Relrow, bool> {
	private:
		float cutoff;
	public:
		explicit Relrow_eq1Matchlen(const float cut) : cutoff(cut) { }
		bool operator()(const Relrow& r) const {
			return r.matchlen == cutoff; 
		}
};

///////////////// end of function objects /////////////////
//
//void trim(list<Relrow> &rel, float cutoff);
// returns the minmal ngidentity
//vector<float> trim(list<Relrow> &rel);
struct Relmin {
	float ngiden, smcov, lgcov;
	int matchlen;
	// these parameter can be passed from the program interface
	bool good() { return (ngiden > 0.5 && smcov > 0.5 && lgcov > 0.7); }
	bool bad() { return !good(); }
};
	
Relmin trim(list<Relrow> &rel);
//void trim(list<Relrow> &rel);
void list2map(const list<Relrow> &rel, multimap<string,string> &mm);
void usage();
//////////////////////////////// Parameters for this program //////////////
struct Subclparam {
	string cltab;  // cluster table
	string alntab;  // alignment table, source of input for relations
	string treetab;  // the tree table: gpcrsp_cltree
	//string covcutoffstr;  // minimal coverage needs to be used in alntab
	string matchlencutoffstr;  // this one seems to be useless
};
// we use two coverage cuttoffs, the smaller and the larger
void breakonecl(connection &pgdb, const string &clidstr, const Subclparam &subparam);

bool UPDATE_DB = false;
set<string> getClusterMembers(const list<Relrow> &rel);
//bool isParent(const string &clid, const Subclparam &subparam, PgDatabase &pgdb);
bool isParent(const string &clid, const Subclparam &subparam, work &TR);

int main(int argc, char* argv[]) {
	// database connection informaiton
	string pghost="hum";
	string pgdbname = "fribi", pguser="kzhou";

	///// program parameters, default settings
	// The can be changed through the program interface
	Subclparam params;
	params.cltab = "gpcrsp_cl";
	params.alntab = "gpcrsp_slfaln";
	params.treetab = "gpcrsp_cltree";
	//params.covcutoffstr = "0.3";
	params.matchlencutoffstr = "100";
	//////////////////////////////////////////

	//string clidstr;
	vector<string> clids;  // in string format
	int i = 1;
	while (i < argc) {
		if (!strcmp(argv[i], "-d")) pgdbname = argv[++i];
		else if (!strcmp(argv[i], "-h")) pghost = argv[++i];
		else if (!strcmp(argv[i], "-u")) pguser = argv[++i];
		else if (!strcmp(argv[i], "--update-db")) UPDATE_DB = true;
		//else if (!strcmp(argv[i], "-of")) resultFile = argv[++i];
		else { // the cluster_id's 
			clids.push_back(argv[i]);
		}
		++i;
	}

	// establish connection to the database
	string pgconnstr = "host=" + pghost + " dbname=" + pgdbname + " user=" + pguser;
	//PgDatabase pgdb(pgconnstr.c_str());
	connection pgdb(pgconnstr.c_str());
   /*
	if (pgdb.ConnectionBad()) {
		cerr << "Connectiion to " << pgconnstr << " failed\n";
		exit(1);
	}
   */
	cerr << "connected to " << pgconnstr << endl;
   work Xaction(pgdb);

	for (int i=0; i<clids.size(); i++) {
		cout << "working on cluster " << clids[i] << " ...\n";
		if (isParent(clids[i], params, Xaction)) {
			cerr << "Cluster is not a leaf, Should not break!\n";
			continue;
		}
		breakonecl(pgdb, clids[i], params);
	}
	return 0;
}

void breakonecl(connection &pgdb, const string &clidstr, const Subclparam &subparam) 
{
	// obtain all relations of input cluster
//	string query = "select query, qcov, target, tcov, matchlen, identity, ngidentity, similarity, score from gpcrsp_slfaln where query in ( select acc from gpcrsp_cl where cluster_id=" + clidstr + ") and target in ( select acc from gpcrsp_cl where cluster_id=" + clidstr + ") and quality >= 0 and (qcov > " + subparam.covcutoffstr + " or tcov > " + subparam.covcutoffstr + " or matchlen > " + subparam.matchlencutoffstr + ") and (qcov>0.12 and tcov>0.12)";

	// No filtering, just select all relations that 
	// produced this cluster, let this program do the filtering
	string query = "select query, qcov, target, tcov, matchlen, identity, ngidentity, similarity, score from "
		+ subparam.alntab
		+ " where query in ( select acc from " + subparam.cltab 
			+ " where cluster_id=" +	clidstr 
		+ ") and target in ( select acc from " 
			+ subparam.cltab + " where cluster_id=" + clidstr 
		+ ") and quality >= 0 and ngidentity >= (select min_ngiden from "
			+ subparam.treetab + " where cluster_id=" + clidstr 
		+ ") and (qcov >= (select min_smcov from "
			+ subparam.treetab + " where cluster_id=" + clidstr
		+ ") and tcov >= (select min_smcov from "
			+ subparam.treetab + " where cluster_id=" + clidstr
		+ ")) and (qcov >= (select min_lgcov from "
			+ subparam.treetab + " where cluster_id=" + clidstr
		+ ") or tcov >= (select min_lgcov from "
			+ subparam.treetab + " where cluster_id=" + clidstr + "))";

   work Xaction(pgdb);
   result qres=Xaction.exec(query);
   /*
	if (!pgdb.ExecTuplesOk(query.c_str())) {
		cerr << pgdb.ErrorMessage() << "Query: " << query << " failed\n";
		exit(1);
	}
   */
	int rows = qres.size();
	int i;
	cerr << "Got " << rows << " rows from " << subparam.alntab << "\n";
	list<Relrow> queryResult;
	for (i=0; i<rows; i++) {
		queryResult.push_back(Relrow(qres, i));
	}
	set<string> parentMembers = getClusterMembers(queryResult);
	Relmin minval = trim(queryResult);
	multimap<string,string> qtpair;
	list2map(queryResult, qtpair);
	hatrees<string> cluster;
	cluster.readFromMap(qtpair);
	vector< set<string> > clarray;
	cluster.clusterArray(clarray);
	// the subclustering may also stop when the whole cluster
	// looks very good
	while (clarray.size() == 1 && minval.bad())
	{
		minval = trim(queryResult);
		list2map(queryResult, qtpair);
		cluster.clear();
		cluster.readFromMap(qtpair);
		cluster.clusterArray(clarray);
	}
	if (clarray.size() == 1) {
		if (parentMembers == getClusterMembers(queryResult)) {
			cerr << "Newly generated cluster same as parent\n"
				<< "Cannot break " << clidstr << endl;
			return;
		}
	}
	// report sime basic information
	cout << "Broken into " << clarray.size() << " subclusters\n"
		<< "with min_ngiden=" << minval.ngiden
		<< " min_smcov=" << minval.smcov
		<< " min_lgcov=" << minval.lgcov << endl;

	// find the largest cluster_id in the database
	query = "select max(cluster_id) from " + subparam.cltab;
   try {
      qres=Xaction.exec(query);
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
   Xaction.commit();
   /*
	if (!pgdb.ExecTuplesOk(query.c_str())) {
		cerr << "faile to execute dbquery: " << query << endl;
		exit(1);
	}
   */

	//int clid_counter = atoi(pgdb.GetValue(0,0)) + 1;
	int clid_counter;
   qres[0][0].to(clid_counter);
   clid_counter++;

	// update the database, cluster table

	string newclstr;
	if (UPDATE_DB) {
		string query_comm = "insert into " + subparam.cltab + "(acc, cluster_id) values(";
		//pgdb.Exec("Begin");
		for (i=0; i<clarray.size(); i++) {
			set<string>::const_iterator it=clarray[i].begin();
			newclstr = itos(i + clid_counter);
			while (it != clarray[i].end()) {
				query = query_comm + "'" + *it + "', " + newclstr + ")";
            try {
               Xaction.exec(query);
            }
            catch (exception &err) {
               //if (!pgdb.ExecCommandOk(query.c_str())) {
					//cerr << pgdb.ErrorMessage() << endl;
					//pgdb.Exec("ROLLBACK");
               cerr << err.what() << endl;
					exit(1);
				}
            Xaction.commit();
				++it;
			}

			// insert into gpcrsp_cltree
			// needs std_, avg_, min_ngiden, min_smcov, min_lgcov
			// compute statistics 
			stddev s;
			float min_ngiden=1, min_smcov=1, min_lgcov=1, relcount=0;

			list<Relrow>::const_iterator li;
			for (li=queryResult.begin(); li != queryResult.end(); li++) {
				if (clarray[i].find(li->query) != clarray[i].end()
						and clarray[i].find(li->target) != clarray[i].end()) 
				{
					s(li->ngidentity);
					if (li->ngidentity < min_ngiden) 
						min_ngiden = li->ngidentity;
					if (li->smcov() < min_smcov)
						min_smcov = li->smcov();
					if (li->lgcov() < min_lgcov)
						min_lgcov = li->lgcov();
					++relcount;
				}
			}
			pair<double,double> avgstd = s.result();

			ostringstream ostr;
			ostr << "insert into " << subparam.treetab
				<< " (cluster_id, parent, min_ngiden, min_smcov, min_lgcov, avg_ngiden, std_ngiden, relcount) values("
				<< newclstr << ", " << clidstr
				<< ", " << min_ngiden
				<< ", " << min_smcov << ", " << min_lgcov
				<< ", " << s.getMean() << ", " << s.getSampleStd()
				<< ", " << relcount << ")";
			string treequery = ostr.str();
         try {
            Xaction.exec(treequery);
         }
         catch (exception &err) {
            cerr << err.what() << endl;
            exit(1);
         }
         Xaction.commit();
         /*
			if (!pgdb.ExecCommandOk(treequery.c_str())) {
				cerr << pgdb.ErrorMessage() << endl
					<< treequery << " Failed!\n";
				pgdb.Exec("ROLLBACK");
				exit(1);
			}
         */
		}
		//pgdb.Exec("COMMIT");
	}
	else {
		string resultFile = "sub" + clidstr + ".cl";
		ofstream OU(resultFile.c_str());
		cluster.showClusterIntId(OU, clid_counter);
		OU.close();
		cerr << "result written to file: " << resultFile << endl;
	}
}

//bool isParent(const string &clid, const Subclparam &subparam, PgDatabase &pgdb) {
bool isParent(const string &clid, const Subclparam &subparam, work &TX) {
	// to prevent from breaking clusters that are already broken
	string query="select * from " + subparam.treetab + " where parent=" + clid;
   result res;
   try {
      res=TX.exec(query);
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
   TX.commit();
   /*
	if (!pgdb.ExecTuplesOk(query.c_str())) {
		cerr << pgdb.ErrorMessage() << endl
			<< query << " Failed!\n";
	}
   */
	if (res.size() > 0) {
		return true;
	}
	return false;
}

set<string> getClusterMembers(const list<Relrow> &rel) {
	set<string> members;
	list<Relrow>::const_iterator it;
	for (it=rel.begin(); it != rel.end(); it++) {
		members.insert(it->query);
		members.insert(it->target);
	}
	return members;
}

void list2map(const list<Relrow> &rel, multimap<string,string> &mm) {
	if (!mm.empty()) mm.clear();
	if (rel.size() < 1) {
		cerr << "List empty in list2map()\n";
		exit(1);
	}

	list<Relrow>::const_iterator it = rel.begin();

	while (it != rel.end()) {
		mm.insert(make_pair(it->query, it->target));
		++it;
	}
}

// remove the edges with the smallest values
// The statistics is a transient state, and should not be 
// saved
// Return a dumy Relrow to hold all the minimums
//vector<float> trim(list<Relrow> &rel) {
//void trim(list<Relrow> &rel) {
Relmin trim(list<Relrow> &rel) {
	if (rel.size() < 2) {
		cerr << "only tow members left! Cannot trim anymore\n";
		exit(1);
	}
	Relmin minval;

	minval.ngiden = min_element(rel.begin(), rel.end(), Relrow_less())->ngidentity;
	cout << "\nminimum: ngidentity = " << minval.ngiden << "; ";
	//it = min_element(rel.begin(), rel.end(), Relrow_lessScore());
	//cerr << "minimum score = " << it->score << endl;
	minval.smcov = min_element(rel.begin(), rel.end(), Relrow_lessSmcov())->smcov();
	cout << "Coverage smaller = " << minval.smcov << ", ";
	minval.lgcov = min_element(rel.begin(), rel.end(), Relrow_lessLgcov())->lgcov();
	cout << "larger = " << minval.lgcov << endl;
	if (minval.lgcov < minval.ngiden) {
		cout << "better trim lgcov\n";
		rel.remove_if(Relrow_eq1Lgcov(minval.lgcov));
		if (minval.smcov < minval.lgcov-0.2) {
			rel.remove_if(Relrow_eq1Smcov(minval.smcov));
		}
	}
	else if (minval.smcov < minval.ngiden-0.2) {
		cout << "better trim smcov\n";
		rel.remove_if(Relrow_eq1Smcov(minval.smcov));
	}
	else {
		cout << "better trim ngiden\n";
		rel.remove_if(Relrow_eq1(minval.ngiden));
	}

	// remove some extremely short matches
	minval.matchlen = min_element(rel.begin(), rel.end(), Relrow_lessMatchlen())->matchlen;
	stddev ss;
	list<Relrow>::iterator it;
	for (it=rel.begin(); it != rel.end(); it++) {
		ss(it->matchlen);
	}
	if ((minval.matchlen < 30 && minval.matchlen < ss.getMean()*0.3) || 
			(minval.matchlen < 100 && minval.matchlen < ss.getMean()*0.18))
	{
		cout << "Removing very short matches: " << minval.matchlen
			<< " avgmatchlen: " << ss.getMean() 
			<< " stdmatchlen: " << ss.getSampleStd() << endl;
		rel.remove_if(Relrow_eq1Matchlen(minval.matchlen));
	}
	//vector<float> clstat;
	//clstat.push_back(minngiden);
	//clstat.push_back(minsmcov);
	//clstat.push_back(minlgcov);
	//return clstat;
	//return minngiden;
	return minval;
}

void usage() {
	cerr << "usage: breakcl [list of cluster_id's]";
	exit(1);
}

