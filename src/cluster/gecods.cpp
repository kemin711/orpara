// file: gecods.cpp
#include "hatrees.h"
//#include <stack>
#include <map>
//#include "/usr/local/pgsql/include/libpq++.h"
//#include "libpq++.h" // no longer avaialble!
#include <pqxx>
#include <string>
#include "loadpgdb.h"
//#include <multimap.h>
#include <iostream>
#include <cstring>

using namespace std;
using namespace pqxx;

/* when runtime error message
 * error while loading shared libraries: libpq++.so.4: cannot open shared object file: No such file or directory
 * either edit /etc/ld.so.conf or define LD_LIBRARY_PATH
 * to include /usr/local/pgsql/lib
 * */

void usage(ostream &os);
enum IdType {STRID, INTID};

int main(int argc, char* argv[]) {
	string infile="fr.slf", pghost = "localhost", pgdb="orpara"; 
	string input_pgtable="self_cov_sel";
	string output_pgtable="same_cds";
	string output_file="same.tab";
	bool fileOutput = true;  // a lot faster for loading
	string inputQuery;
	bool fileInput = false;
	bool coverage_cut = false;
	string cut = "0.8";   // coverage cut, used for sql, sting is better
	string key_type = "string";  // string or int
	IdType idtype = INTID; // 0 for string as id, 1 for integer as id
	int clid = 1;
	
	if (argc < 2) usage(cerr);
	int i=1;
	while (i<argc) {
		if (!strcmp(argv[i], "-f")) {
			infile=argv[++i];
			fileInput = true;
			fileOutput = true;
		}
		else if (!strcmp(argv[i], "-h")) pghost=argv[++i];
		else if (!strcmp(argv[i], "-d")) pgdb=argv[++i];
		else if (!strcmp(argv[i], "-i")) input_pgtable=argv[++i];
		else if (!strcmp(argv[i], "--idtype")) {
			++i;
			if (!strcmp(argv[i],"int")) {
				idtype = INTID;
			}
			else if (!strcmp(argv[i], "str")) idtype = STRID;
			else {
				cerr << "idtype must be either int or str\n";
				exit(1);
			}
		}
		else if (!strcmp(argv[i], "--clid")) {
			clid = atoi(argv[++i]);
			idtype = INTID;
		}
		else if (!strcmp(argv[i], "-o")) { 
			fileOutput = false;
			output_pgtable=argv[++i];
		}
		else if (!strcmp(argv[i], "-of")) {
			if (++i >= argc) {
				cerr << "you must supply a fileName after -of option\n";
				exit(1);
			}
			output_file=argv[i];
		}
		else if (!strcmp(argv[i], "-c")) {
			cut = argv[++i];
			coverage_cut = true;
		}
		else if (!strcmp(argv[i], "-q")) inputQuery = argv[++i];
		else if (!strcmp(argv[i], "-k")) key_type = argv[++i];
		else { 
			cerr << argv[i] << " is not a legal option\n";
			usage(cerr); 
		}
		i++;
	}

	hatrees<string> clusters;
	ofstream OUT;
	if (fileOutput) {
		OUT.open(output_file.c_str());
		if (OUT.fail()) { exit(1); }
	}
	if (fileInput) {
		clusters.readFromFile(infile);
		OUT << "RESULTS: " << clusters.getNodeCount() << "\n";
		clusters.showStore(OUT);
		cout << "Now it is sorted according to cluster id\n";
		clusters.showCluster(OUT);
		clusters.showClusterByLine(cout);
		return 0;
	}

	////// working from a database /////////////////////
	string connInfo = "host=" + pghost + " dbname=" + pgdb;
	//PgDatabase orthodb(connInfo.c_str());
	connection orthodb(connInfo.c_str());
   /*
	if (orthodb.ConnectionBad()) {
		cerr << "Connectiion to " << connInfo << " failed\n";
		exit(1);
	}
   */
	cerr << "Connected to " << connInfo << endl;

	if (inputQuery.empty()) { // will make query
		if (coverage_cut) {
			inputQuery = "select query, target from ";
			inputQuery.append(input_pgtable).append(" where qcov>");
			inputQuery.append(cut).append(" or tcov > ").append(cut);
		}
		else inputQuery = "select query,target from " + input_pgtable;
	}
   work Xaction(orthodb);
   result qres=Xaction.exec(inputQuery);
   /*
	if(!orthodb.ExecTuplesOk(inputQuery.c_str())) {
		cerr << "Query: " << inputQuery << " failed\n";
		exit(1);
	}
   */

	hatrees<int> int_clusters;

	if (key_type == "string") {
		//clusters.readFromDB(orthodb);
		clusters.readFromDB(qres);
		cerr << "Key type: string; read from DB done.\n";
		clusters.showClusterByLine(cout);
		cout << clusters.getNodeCount() << " unique members\n";
		if (fileOutput)  {
			if (idtype == INTID) clusters.showClusterIntId(OUT, clid);
			else clusters.showCluster(OUT); 
		}
		else clusters.loadTable(orthodb, output_pgtable);
	}
	else {
		//int_clusters.readFromDB(orthodb);
		int_clusters.readFromDB(qres);
		int_clusters.showClusterByLine(cout);
		if (fileOutput)  int_clusters.showCluster(OUT); 
		else int_clusters.loadTable(orthodb, output_pgtable);
	}

	return 0;
}

void usage(ostream &os) {
	os << "Usage: gecods -f infile \n"
		<< "\t-d pgdb           default=orpara\n"
		<< "\t-h pghost         default=localhost\n"
		<< "\t-i input_pgtable  default=self_cov_sel\n"
		<< "\t-o output_pgtable default=same_cds\n"
		<< "\t-of output_file   will only output to file\n"
		<< "\t-c 0.1--1.0       default=0.8 optional\n"
		<< "\t-q [query]        default=select query,target from input_table\n"
		<< "\t-k key_type       int or string  default=string\n"
		<< "\t--idtype str or int if --idtype set to int then needs clid\n"
		<< "\t--clid integer    provide the initial cluster id default=1\n"
		<< "\t\toptional\n";
	os << "file input or db input are exclusive\n";
	exit(1);
}
