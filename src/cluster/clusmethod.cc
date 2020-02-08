// file: clusmethod.cpp
//
#include "clusmethod.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

//char clusmethod::connInfo[15] = "dbname=ortho";

/*
clusmethod::clusmethod(const char *query, const char *connInfo) 
	: orthodb(connInfo) 
{
	if (orthodb.ConnectionBad()) {
		cerr << "Connecting to " << connInfo << " db failed.\n";
		exit(1);
	}
	cerr << "Connected to " << connInfo << "\n";
	//char statement[] = "SELECT query, target FROM cds_dna_homol";
	//char statement2[] = "SELECT target, query FROM cds_dna_homol ORDER BY target";
	if (!orthodb.ExecTuplesOk(query)) {
		cerr << "SQL: " << query << " failed\n";
		exit(1);
	}
	int i, first, second;
	for (i=0; i<orthodb.Tuples(); i++) {
		first = atoi(orthodb.GetValue(i,0));
		second = atoi(orthodb.GetValue(i,1));
		rel.insert( pair<int, int>(first, second) );
:q
		revrel.insert( pair<int, int>( second, first) );
		//rel.insert( pair<int, int>( atoi(orthodb.GetValue(i,0)),  atoi(orthodb.GetValue(i,1))));
		//revrel.insert( pair<int, int>( atoi(orthodb.GetValue(i,1)),  atoi(orthodb.GetValue(i,0))));
	}
}
*/

/* simple interface, query result must contain (query,target) */

//void clusmethod::dbinput(PgDatabase &orthodb) 
void clusmethod::dbinput(result &qres) 
{
	int first, second;
   result::size_type i;
	for (i=0; i<qres.size(); i++) {
		//first = atoi(qres[i][0]);
      qres[i][0].to(first);
		//second = atoi(orthodb.GetValue(i,1));
		//second = atoi(qres[i][1]);
		qres[i][1].to(second);
		rel.insert( pair<int, int>(first, second) );
		revrel.insert( pair<int, int>( second, first) );
	}
}

void clusmethod::fileinput(const char file[]) {
	// not database connection needed
	// to be implemented
	// a problem with PgDatabase! it does not allow you to 
	// costruct a empty database
	// file with
	// first_number second_number separator=\t
	ifstream IN(file);
	if (IN.fail()) {
		cerr << "File: " << file << " open FAILED\n";
		exit(1);
	}
	int first, second;
	IN >> first >> second;
	while (!IN.eof())  {
		rel.insert(make_pair(first, second));
		revrel.insert(make_pair(second, first));
		IN >> first >> second;
	}
}

vector<cluster>& clusmethod::findCluster() {
	int seedKey;
	set<int> orthos, newOrthos;

	while (rel.size() > 0) {
		allcl.push_back(cluster(rel));
		newOrthos = processTargets(allcl.back().set2);
		while (newOrthos.size() > 0) {
			allcl.back().addTargets(newOrthos);
			orthos = newOrthos;
			newOrthos = processTargets(orthos);
		}
	}
	return allcl;
}

set<int> clusmethod::processTargets(const set<int> &tset) {
	set<int> newTargets;  // to find some new target for our cluster
	set<int>::iterator tseti = tset.begin();

	while (tseti != tset.end()) {
		int targetKey = *tseti;
		set<int> tmpSeedSet = grab(targetKey, revrel);
		set<int>::iterator sdSeti = tmpSeedSet.begin();
		while (sdSeti != tmpSeedSet.end()) {
			int nextSeed = *sdSeti;
			if (rel.find(nextSeed) != rel.end()) {
				allcl.back().addSeed(nextSeed);
				set<int> moreTargets = grab(nextSeed, rel);
				set_difference(moreTargets.begin(), moreTargets.end(), 
						allcl.back().set2.begin(), allcl.back().set2.end(),
						insert_iterator<set<int> >(newTargets, newTargets.begin()) );
			}
			sdSeti++;
		}
		tseti++;
	}
	return newTargets;
}

void clusmethod::printMap()
{
	multimap<int, int>::iterator i = rel.begin();
	while (i != rel.end()) {
		cout << i->first << " " << i->second << endl;
		i++;
	}
	cout << "\n============ reverse map ============\n";
	i = revrel.begin();
	while (i != revrel.end()) {
		cout << i->first << " " << i->second << endl;
		i++;
	}
}

/**
 * helper function used only by processTargets
 * It is declared as friend but is giving problems.
 */
set<int> grab(int k, multimap<int, int> &mm) {
	typedef multimap<int, int>::iterator mmit;
	set<int> tmp;

	pair<mmit, mmit> bp = mm.equal_range(k);
	mmit p = bp.first;
	if (p != mm.end()) { // found key k
		while (p != bp.second) {
			tmp.insert(p->second);
			p++;
		}
		mm.erase(k);
	}
	return tmp;
}

void clusmethod::displayClusters(ostream &os) {
	os << "seeds => targets\n--------------------\n";
	for (int i=0; i<allcl.size(); i++) {
		os << allcl[i] << endl;
	}
}

/* col is the column name of the cds table.  default is 
 * n_cluster_id; cluster_id is the starting number and will be 
 * incremented for each cluster added*/
bool clusmethod::loadTable(connection &orthodb, string update_tab, 
		string update_col, string clutab, string clus_col, int &cluster_id) 
{
	set<int> allMembers;
	set<int>::iterator si;
	string arrStr, query; 
	//string col_cluster = "cds_set";  // col into which to insert 
	string queryStart = "INSERT INTO " +  clutab;
	queryStart.append(clus_col + " values(");

	string cdsUpdateHead = "UPDATE " + update_tab + " set " + update_col + "=";

	string cdsUpdate;
	const int CLUSTER_SIZE_LIMIT=1000;  // anything larger than 1000 is 
	// not inserted into the insert_table
   work Xaction(orthodb, "forinsertion");

	for (int i=0; i<allcl.size(); i++) {
		cerr << "working on " << i << "th cluster, with id "
			<< cluster_id << endl;
		ostringstream conv;
		conv << cluster_id << ", '";
		if (allcl[i].size() > CLUSTER_SIZE_LIMIT) {
			conv << "Large, with " << allcl[i].size() << " members";
		}
		else {
			conv << allcl[i].toText();
		}
		conv << "')";
		query = queryStart + conv.str();
      try {
         Xaction.exec(query);
      }
      catch (exception &err) {
			cerr << "Failed to insert into " << err.what() << endl;
			return false;
		}
      Xaction.commit();
		/* the insertion process is very slow in the new version of the
		 * database; It involves a lot disk access. Better use dumpTable
		 * and load the table with \copy command
		 */

		allMembers = allcl[i].getElements();
		cerr << "\t" << allMembers.size() << " members\n";
		si = allMembers.begin();
		while (si != allMembers.end()) {
			ostringstream converter;
			converter << cluster_id << " where id=" << *si;
			cdsUpdate = cdsUpdateHead + converter.str();
         try {
            Xaction.exec(cdsUpdate);
         }
         catch (exception &err) {
				cerr << "Failed to insert into " << update_tab 
					<< " with the values: " << err.what() << endl;
				return false;
			}
			si++;
		}
		cluster_id++;
	}
	return true;
}
void clusmethod::dumpTable(const char ins[], const char upd[], int &cluster_id)
{
	set<int> allMembers;
	set<int>::iterator si;

	ofstream INS(ins);
	ofstream UPD(upd);
	if (INS.fail() || UPD.fail()) {
		cerr << "Failed to open " << ins << " or " << upd << " file\n";
		exit(1);
	}

	for (int i=0; i<allcl.size(); i++) {
		cerr << "working on " << cluster_id << " cluster\n";
		if (allcl[i].size()>1000) {
			INS << cluster_id << '\t' << "Large cluster: " << allcl[i].size()
				<< " members\n";
		}
		else {
			INS << cluster_id << '\t' << allcl[i].toText(' ') << endl;
		}

		allMembers = allcl[i].getElements();
		cerr << "\t" << allMembers.size() << " members\n";
		cout << cluster_id << '\t' << allMembers.size() << endl;
		si = allMembers.begin();
		while (si != allMembers.end()) {
			UPD << *si << '\t' << cluster_id << endl;
			si++;
		}
		cluster_id++;
	}
	INS.close();
	UPD.close();
}
