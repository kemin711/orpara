/* this program use the input from the output of clustalw
 * and do subcluster.  uset the *.score file 
 * This is not ncessary because I have invented a much better way
 * of doing the same job using SQL and results from tablan
 * */

//#include "forest.h"
#include "hatrees.h"
#include <map>
#include <set>
#include <algorithm>
#include <unistd.h>
#include <dirent.h>
#include "scorepair.h"
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <cstring>

using namespace std;


enum quality { good, bad, modify };

void read(ifstream &in, multimap<int, int> &m, int cut);
quality qual(ifstream &in, int cut);
void usage();

int main(int argc, char* argv[]) {
	bool loaddb = true;
	int i = 1;
	string infile, indir="."; //indir: input directory relative to current
	int scoreCut = 23;
	bool doall = true;
	while (i < argc) {
		if (!strcmp(argv[i], "-c")) scoreCut = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-d")) indir = argv[++i];
		else if (!strcmp(argv[i], "-h")) usage();
		else infile = argv[i];
		i++;
	}
	int clusterIdCnt;
	char conninfo[] = "host=lcc2 dbname=ortho user=o_admin";
	//PgDatabase pgdb(conninfo);
	connection pgdb(conninfo);
   work Xaction(pgdb, "forqueryPGDB");
	if (loaddb) {    
      /*
		if (pgdb.ConnectionBad()) {
			cerr << conninfo << " failed\n";
			exit(1);
		}
      */
		//pgdb.ExecTuplesOk("select max(cluster_id) from porthopara");
		result qres=Xaction.exec("select max(cluster_id) from porthopara");
      Xaction.commit();
		//clusterIdCnt = atoi(pgdb.GetValue(0,0));
		qres[0][0].to(clusterIdCnt);
		cout << "last id used " << clusterIdCnt << endl;
	}

	if (!infile.empty()) {
		ifstream IN(infile.c_str());
		if (IN.fail()) {
			cerr << infile << " opening failed\n";
			exit(1);
		}
		scorepair subcl;
		subcl.read(IN, 23);
		int q = subcl.quality();
		cout << "quality is " << q << endl;
		IN.close();
	}
	else {  // doing all *.score files
		//char *cwd = new char[501];
		//getcwd(cwd, 500);
		cout << " working on directory: " << indir << endl;
		int oo;
		if (indir != ".") { oo = chdir(indir.c_str()); }
		if (oo == -1) { 
			cerr << "chdir " << indir << " failed\n"; 
			exit(1);
		}
		DIR *curdir = opendir(".");
		dirent *afile = readdir(curdir);
		//quality qqq;
		int qqq;
		int goodCnt=0, badCnt=0, modifyCnt=0;
		string clusterid;
		while (afile != NULL) {
			char *p = strstr(afile->d_name, ".score");
			if (p && *(p+6) == '\0') {
				clusterid = string(afile->d_name, p-afile->d_name);

				ifstream IN(afile->d_name);

				scorepair sp;
				sp.read(IN, scoreCut);
				qqq = sp.quality();

				if (loaddb) {
					cout << " load result into database ...\n";
					string query;
					if (qqq == 1) { // good
						goodCnt++;
						query = "update porthopara set quality='3' where cluster_id=";
						query += clusterid;
                  Xaction.exec(query);
                  /*
						if (!pgdb.ExecCommandOk(query.c_str())) {
							cerr << query << " failed\n";
							exit(1);
						}
                  */
					}
					else if (qqq = 9) {
						deleteCluster(pgdb, clusterid);
						badCnt++;
					}
					else {
						modifyCnt++;
						cout << " modifying cluster " << clusterid << endl;
						sp.addNewCluster(pgdb, clusterid, clusterIdCnt);
					}
				}
				else {
					if (qqq == 1) {
						cout << clusterid << " good\n";
						goodCnt++;
					}
					else if (qqq == 9) {
						cout << clusterid << " bad\n";
						badCnt++;
					}
					else {
						cout << clusterid << " modified\n";
						modifyCnt++;
						sp.showcluster(cout);
						set<string> bdmb = sp.badMembers();
						cout << "bad members: \n";
						copy(bdmb.begin(), bdmb.end(), ostream_iterator<string>(cout, " "));
					}
				}

				cout << endl;
				IN.close();
			}
			afile = readdir(curdir);
		}
		closedir(curdir);
		//delete[] cwd;
		cout << "good: " << goodCnt << " bad: " << badCnt << " modified: "
			<< modifyCnt << endl;
	}

	return 0;
}

void usage() {
	cerr << "usage: subcluster, default behavior will working on "
		<< "\nthe whole directory of clustal results\n";
	cerr << "usage: subcluster -c [score_cut default 23] infile\n"
		<< " -h for help\n";
	exit(1);
}

/** read the score file, discarding any pair with 
 * score < cut */
void read(ifstream &in, multimap<int, int> &m, int cut) {
	static const int LN = 100;
	char ln[LN+1];
	string name, dumy;
	int len, snum;
	vector<string> seqnames;
	vector<int> seqlens;

	in.getline(ln, LN);
	while (strncmp(ln, "Sequence 1", 10)) in.getline(ln, LN);
	//cout << ln << endl;
	len = strlen(ln);
	//in.seekg(-(strlen(ln)+1), ios::cur);  // this one just won't work!
	in.seekg(-len, ios::cur);             // this works! Why?
	in >> ln >> dumy >> name >> len >> dumy;
	//cout << "name: " << name << " len:" << len << endl;
	seqnames.push_back(name);
	seqlens.push_back(len);   // information not used yet
	in >> ln;
	while (!strcmp(ln, "Sequence")) {   // read name length lines
		in >> dumy >> name >> len >> dumy;
		//cout << name << " " << len << endl;
		seqnames.push_back(name);
		seqlens.push_back(len);
		in >> ln;
	}
	char *p, *pp;
	int f1, f2, score;
	m.clear();

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
		if (score > cut) {
			m.insert(pair<int, int>(f1, f2));
		}
		else {
			cout << f1 << " " << f2 << " discarded " << score << endl;
		}
		in.getline(ln, LN);
	}
}

/* will do everything on the input score file */
quality qual(ifstream &in, int cut) {
	static const int LN = 100;
	char ln[LN+1];
	string name, dumy;
	int len, snum;
	vector<string> names;
	vector<int> seqlens;
	multimap<string, string> m;
	//names.clear();

	in.getline(ln, LN);
	while (strncmp(ln, "Sequence 1", 10)) in.getline(ln, LN);
	len = strlen(ln);
	in.seekg(-len, ios::cur);             // this works! Why?
	in >> ln >> dumy >> name >> len >> dumy;
	names.push_back(name);
	seqlens.push_back(len);   // information not used yet
	in >> ln;
	while (!strcmp(ln, "Sequence")) {
		in >> dumy >> name >> len >> dumy;
		//cout << name << " " << len << endl; // debug
		names.push_back(name);
		seqlens.push_back(len);
		in >> ln;
	}
	//cout << names.size() <<  " sequences\n";
	char *p, *pp;
	int f1, f2, score;
	//m.clear();  // done inside this function
	//set<string> discard;  // discarded sequences names

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
		if (score > cut) {
			m.insert(pair<string, string>(names[f1-1], names[f2-1]));
		}
		in.getline(ln, LN);
	}
	if (m.empty()) { cout << "cluster bad\n"; return bad; }
	else {
		hatrees<string> subclstr;
		subclstr.readFromMap(m);
		if (names.size() == subclstr.getNodeCount()) {
			cout << "cluster good\n";
			return good;
		}
		else {
			cout << "total " << names.size() << " sequences\n";
			set<string> retained = subclstr.keyset();
			set<string> discard;
			sort(names.begin(), names.end());
			set_difference(names.begin(), names.end(), retained.begin(), retained.end(), inserter(discard, discard.begin()));
			cout << "discarded ones: ";
			copy(discard.begin(), discard.end(), ostream_iterator<string>(cout, " "));
			cout << endl;
			//copy(retained.begin(), retained.end(), ostream_iterator<string>(cout, " "));
			cout << endl;
			subclstr.showClusterByLine(cout);
			return modify;
		}
	}
}
