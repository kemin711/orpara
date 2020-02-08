#include "scorepair.h"
#include <sstream>
#include <iostream>
#include <cstring>
#include <algorithm>


void scorepair::read(istream &in, int cut) {
	clear();
	static const int LN = 100;
	char ln[LN+1];
	string name, dumy;
	int len, snum;
	//vector<string> names;
	//vector<int> seqlens;
	//multimap<string, string> m;

	in.getline(ln, LN);
	while (strncmp(ln, "Sequence 1", 10)) in.getline(ln, LN);
	len = strlen(ln);
	in.seekg(-len, ios::cur);             // this works! Why?
	in >> ln >> dumy >> name >> len >> dumy;
	names.push_back(name);
	seqlens.push_back(len);   // information not used yet
	in >> ln;
	while (!strcmp(ln, "Sequence")) {  // read sequence names
		in >> dumy >> name >> len >> dumy;
		names.push_back(name);
		seqlens.push_back(len);
		in >> ln;
	}
	char *p, *pp;
	int f1, f2, score;

	in.getline(ln, LN);
	while (strncmp(ln, "Sequences", 9)) in.getline(ln, LN);
	while (!strncmp(ln, "Sequences", 9)) {
		p = ln + 11;
		pp = p + 1;
		while (*pp != ':') pp++;
		*pp = '\0';
		f1 = atoi(p);
		p = pp + 1;
		pp = p + 1;
		while (*pp != ')') pp++;
		*pp = '\0';
		f2 = atoi(p);
		p = pp + 1;
		while (!isdigit(*p)) p++;
		score = atoi(p);
		//cout << names[f1-1] << " " << names[f2-1] << endl;
		if (score > cut) {
			rel.insert(pair<string, string>(names[f1-1], names[f2-1]));
		}
		in.getline(ln, LN);
	}
}

int scorepair::quality() {
	if (rel.empty()) { 
		//cout << "cluster bad\n"; 
		return 9; 
	}

	//forest<string> subclstr;
	subclstr.readFromMap(rel);
	if (names.size() == subclstr.getNodeCount()) {
		//cout << "cluster good\n";
		return 1;
	}
	
	return 5;  // modified
}

void scorepair::showcluster(ostream &ou) {
	subclstr.showClusterByLine(ou);
}

vector<set<string> > scorepair::getClusters() {
	vector<set<string> > vs;
	subclstr.clusterArray(vs);
	return vs;
}

void scorepair::clusterArray(vector<set<string> > &clarr) {
	subclstr.clusterArray(clarr);
}

set<string> scorepair::badMembers() const {
	//cout << "total " << names.size() << " sequences\n";
	set<string> retained = subclstr.keyset();
	set<string> discard;
	set<string> all(names.begin(), names.end());

	//sort(names.begin(), names.end());
	//set_difference(names.begin(), names.end(), retained.begin(), retained.end(), inserter(discard, discard.begin()));
	set_difference(all.begin(), all.end(), retained.begin(), retained.end(), inserter(discard, discard.begin()));
	//cout << "discarded ones: ";
	//copy(discard.begin(), discard.end(), ostream_iterator<string>(cout, " "));
	//cout << endl;
	return discard;
}

void deleteCluster(connection &pgdb, const string &clusterid) {
	string cmd = "update cds set pcluster=null where pcluster=";
	cmd += clusterid;
   work Xaction(pgdb);
   try {
      Xaction.exec(cmd);
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
   Xaction.commit();
   /*
	if (!pgdb.ExecCommandOk(cmd.c_str())) {
		cerr << cmd << " failed\n";
		exit(1);
	}
   */

	cmd = "delete from porthopara where cluster_id=";
	cmd += clusterid;
   try {
      Xaction.exec(cmd);
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
   /*
	if (!pgdb.ExecCommandOk(cmd.c_str())) {
		cerr << cmd << " failed\n";
		exit(1);
	}
   */
}

/* clusterIdCnt is the sequence reading or the max cluster_id from
 * porthopara; must call quality before calling this one */
void scorepair::addNewCluster(connection &pgdb, const string &clusterid,
		int &clusterIdCnt) {
	string query = "select prt, cds_id from cds where pcluster=";
	query += clusterid;
   work Xaction(pgdb);
   result qres;
   try {
      qres=Xaction.exec(query);
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
   Xaction.commit();

   /*
	if (!pgdb.ExecTuplesOk(query.c_str()) ) { 
		cerr << query << " failed\n";
		exit(1);
	}
   */
	map<string, string> prt2id;
	map<string, string>::const_iterator mi;

	int i;
   string tmp1,tmp2;
	//for (i=0; i<pgdb.Tuples(); i++) {
	for (i=0; i<qres.size(); i++) {
		//prt2id[string(pgdb.GetValue(i,0))] = string(pgdb.GetValue(i, 1));
      qres[i][0].to(tmp1);
      qres[i][1].to(tmp2);
		prt2id[tmp1] = tmp2;
	}
	if (prt2id.empty()) {
		cerr << "empty map of prt to cds_id\n";
		exit(1);
	}
	/* once the map is loaded the cluster can be deleted */
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
      try {
         Xaction.exec(query);
      }
      catch (exception &err) {
         cerr << err.what() << endl;
         exit(1);
      }
      Xaction.commit();
      /*
		if (!pgdb.ExecCommandOk(query.c_str())) {
			cerr << query << "failed\n";
			exit(1);
		}
      */

		for (int j=0; j<cds_ids.size(); j++) {
			ostringstream oss1;
			oss1 << "update cds set pcluster=" << clusterIdCnt << " where cds_id="
				<< cds_ids[j];
			query = oss1.str();
         try {
            Xaction.exec(query);
         }
         catch (exception &err) {
            cerr << err.what() << endl;
            exit(1);
         }
         /*
			if (!pgdb.ExecCommandOk(query.c_str())) {
				cerr << query << "failed\n";
				exit(1);
			}
         */
		}
	}
	cout << clusterid << " replaced with " << clarr.size() << " new clusters\n";
}

