#include "loadpgdb.h"
#include <sstream>
#include <iostream>

using namespace std;

//void deleteCluster(PgDatabase &pgdb, const string &clusterid) {
void deleteCluster(pqxx::connection &pgdb, const string &clusterid) {
	string cmd = "update cds set pcluster=null where pcluster=";
	cmd += clusterid;
   pqxx::work Xaction(pgdb);
   Xaction.exec(cmd);
   /*
	if (!pgdb.ExecCommandOk(cmd.c_str())) {
		cerr << cmd << " failed\n";
		exit(1);
	}
   */
   Xaction.commit();

	cmd = "delete from porthopara where cluster_id=";
	cmd += clusterid;
   Xaction.exec(cmd);
   /*
	if (!pgdb.ExecCommandOk(cmd.c_str())) {
		cerr << cmd << " failed\n";
		exit(1);
	}
   */
   Xaction.commit();
}

//void loadClusterIntoDB(PgDatabase &pgdb, 
void loadClusterIntoDB(pqxx::connection &pgdb, 
		const multimap<int, int> &clusters, int &clusterIdCnt) {
	// update cds table and insert into porthopara table
	const string updatetab = "cds";
	const string inserttab = "porthopara";
	string updateQueryBase = "update " + updatetab + " set pcluster=";
	string insertQueryBase = "insert into " + inserttab + " (cluster_id) values(";
	string updateQuery, insertQuery;
   pqxx::work Xaction(pgdb);

	multimap<int, int>::const_iterator mmi = clusters.begin();
	while (mmi != clusters.end()) {
		int cid = mmi->first;
		//insert into porthopara (cluster_id) values(clusterIdCnt);
		ostringstream ssi;  // for inserting cluster_id into porthopara
		ostringstream oou; // for updating cds_set
		ssi << insertQueryBase << clusterIdCnt << ")";

      Xaction.exec(ssi.str());
      Xaction.commit();
      /*
		if (!pgdb.ExecCommandOk(ssi.str().c_str())) {
			cerr << pgdb.ErrorMessage();
			exit(1);
		}
      */
		oou << "update " << inserttab << " set cds_set='";

		while (mmi != clusters.end() && cid == mmi->first) {
			ostringstream ss;
			ss << updateQueryBase << clusterIdCnt << " where cds_id="
				<< mmi->second;
			updateQuery = ss.str();
         Xaction.exec(updateQuery);
         /*
			if (!pgdb.ExecCommandOk(updateQuery.c_str())) {
				cerr << updateQuery << " failed\n";
				exit(1);
			}
         */
         Xaction.commit();
			if (mmi == clusters.begin()) oou << mmi->second;
			else oou << '\t' << mmi->second;
			mmi++;
		}
		oou << "' where cluster_id=" << clusterIdCnt;
      Xaction.exec(oou.str());
      Xaction.commit();
      /*
		if (!pgdb.ExecCommandOk(oou.str().c_str())) {
			cerr << "cds_set update failed at cluster_id=" << clusterIdCnt 
				<< endl;
			exit(1);
		}
      */
		clusterIdCnt++;
	}
}


/* clusterIdCnt is the sequence reading or the max cluster_id from
 * porthopara; must call quality before calling this one */
/*  originally member functions from scorepair.h 
void addNewCluster(PgDatabase &pgdb, const string &clusterid,
		int &clusterIdCnt) {
	string query = "select prt, cds_id from cds where pcluster=";
	query += clusterid;
	if (!pgdb.ExecTuplesOk(query.c_str()) ) { 
		cerr << query << " failed\n";
		exit(1);
	}
	map<string, string> prt2id;
	map<string, string>::const_iterator mi;

	int i;
	for (i=0; i<pgdb.Tuples(); i++) {
		prt2id[string(pgdb.GetValue(i,0))] = string(pgdb.GetValue(i, 1));
	}
	if (prt2id.empty()) {
		cerr << "empty map of prt to cds_id\n";
		exit(1);
	}
	// once the map is loaded the cluster can be deleted 
//
	deleteCluster(pgdb, clusterid);

	vector<set<string> > clarr; 
	clusterArray(clarr);
	set<string>::const_iterator si;
	for (i=0; i<clarr.size(); i++) {
		si = clarr[i].begin();

		ostringstream oss;
		vector<string> cds_ids;
		oss << "insert into porthopara (cluster_id, cds_set) values("
			<< ++clusterIdCnt << ", '";

		mi = prt2id.find(*si);
		if (mi == prt2id.end()) { 
			cerr << "Protein: " << *si << " not found in cluster: " << clusterid
				<< endl;
			exit(1);
		}
		oss << mi->second;
		cds_ids.push_back(mi->second);
		si++;
		while (si != clarr[i].end()) {
			mi = prt2id.find(*si);
			if (mi == prt2id.end()) {   // safer code but ugly
				cerr << "Protein: " << *si << " not found in cluster: " << clusterid
					<< endl;
				exit(1);
			}
			oss << '\t' << mi->second;
			cds_ids.push_back(mi->second);
			si++;
		}
		oss << "')";
		query = oss.str();
		if (!pgdb.ExecCommandOk(query.c_str())) {
			cerr << query << "failed\n";
			exit(1);
		}

		for (int j=0; j<cds_ids.size(); j++) {
			ostringstream oss1;
			oss1 << "update cds set pcluster=" << clusterIdCnt << " where cds_id="
				<< cds_ids[j];
			query = oss1.str();
			if (!pgdb.ExecCommandOk(query.c_str())) {
				cerr << query << "failed\n";
				exit(1);
			}
		}
	}
	cout << clusterid << " replaced with " << clarr.size() << " new clusters\n";
}
*/
