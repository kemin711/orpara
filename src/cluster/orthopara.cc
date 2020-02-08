#include "clusmethod.h"
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>

// file: orthopara.cc
/* this is the main program that used the clusmethod 
 *
 * to to clustering regardless of coverage use the following command:
 * orthopara -q "select query, target from final_bln"  -i 1 
 *    --cc "(cluster_id, cds_set)" --ct northopara1 -u cds1 
 *    -c n_cluster_id
 *
 * to do clustering in separate coverage ranges
 * orthopara -q "select query,target from final_bln where cov>0.1" 
 *    --cut 0.1 i 1
 *
 * */
using namespace std;

void usage(map<string, string> &param);

int main(int argc, char* argv[])
{
	int i=1;

	string host = "localhost";
	string dbname = "orpara";
	string user = "kzhou";
	//////// input table must contain query,target //////////////
	string input_tab = "frgoodpair";
	string query = "select query, target from " + input_tab;
	// query from command line overwrites default query

	/* output section,  updating the cds table */
	string update_tab = "cds";  // updates the cds table: col 
	string update_col = "p_cluster";  // column of the update_tab

	/////////// insert new clusters into the cluster table ////////
	string cluster_tab = "pcluster";
	string cluster_col = "(id, members)";  // columns in the pcluster table
	int cluster_id = 1;  // the first id number for cluster, counter
	bool show = false; // load table, if true, display result to cout
	bool FILEOUT = true; // default shold be file output safer
	string infile;

	map<string, string> progparam;
	progparam["host"]="localhost";
	progparam["dbname"]="orpara";
	progparam["user"]="kzhou";
	progparam["input_tab"]="frgoodpair";
	progparam["update_tab"]="cds";
	progparam["update_col"]="p_cluster";
	progparam["cluster_tab"]="pcluster";
	progparam["cluster_col"]="(id, members)";
	progparam["cluster_id"] = "1";
	progparam["query"]=query;

	if (argc == 1) { 
		usage(progparam); return 1;
	}

	while (i<argc) {
		if (!strcmp(argv[i], "-h")) host = argv[++i];
		else if (!strcmp(argv[i], "-d")) dbname = argv[++i];
		else if (!strcmp(argv[i], "-u")) user = argv[++i];
		else if (!strcmp(argv[i], "-i")) cluster_id = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-it")) {
			input_tab = argv[++i];
			query = string("select query, target from ") + input_tab;
		}
		//else if (!strcmp(argv[i], "-t")) intab = argv[++i];
		else if (!strcmp(argv[i], "-t")) update_tab = argv[++i];
		else if (!strcmp(argv[i], "-c")) update_col = argv[++i];
		else if (!strcmp(argv[i], "-ct")) cluster_tab = argv[++i];
		else if (!strcmp(argv[i], "-cc")) cluster_col = argv[++i];
		//else if (!strcmp(argv[i], "--cut")) cov_cut = argv[++i];
		else if (!strcmp(argv[i], "-show")) show = true;
		else if (!strcmp(argv[i], "-q")) query = argv[++i];
		else if (!strcmp(argv[i], "-f")) infile = argv[++i];
		else if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "--help")) {
			usage(progparam); return 1;
		}
		else if (!strcmp(argv[i], "-fo")) FILEOUT = true;
		else if (!strcmp(argv[i], "-FO")) FILEOUT = false;
		else {
			cerr << argv[i] << " not a legal option\n";
			exit(1);
		}
		i++;
	}

	cerr << "The work summary:\n"
		<< "Starting cluster_id " << cluster_id << endl
		<< "Update " << update_col << " in " << update_tab << endl
		<< "Insert into " << cluster_tab << " Column " <<  cluster_col 
		<< endl
		<< " show=" << show << endl
		<< "query: " <<  query << endl;

	clusmethod theMethod;
	if (infile.empty()) {
		string connectionString = "host=" + host + " dbname="+dbname + " user=" + user;
		//PgDatabase pgdb(connectionString.c_str());
      pqxx::connection pgdb(connectionString.c_str());
      /*
		if (pgdb.ConnectionBad()) {
			 cerr << "Connecting to " << connectionString << " db failed.\n";
			 exit(1);
		}
      */
		cerr << "Connected to " << connectionString << "\n";
      work Xaction(pgdb, "transactionObject");
      result qresult;
      try {
         qresult=Xaction.exec(query);
         /*
         if (!pgdb.ExecTuplesOk(query.c_str())) {
            cerr << "SQL: " << query << " failed\n";
            exit(1);
         }
         */
      }
      catch (exception &err) {
         cerr << err.what() << endl;
         exit(1);
      }
      Xaction.commit();
		theMethod.dbinput(qresult);
		theMethod.findCluster();

		if (FILEOUT) {
			theMethod.dumpTable("pcluster.tab", "cdsupdate.tab", cluster_id);
		}
		else {
			if (!theMethod.loadTable(pgdb, update_tab, update_col, cluster_tab, cluster_col, cluster_id)) 
			{
				cerr << "Failed to load cluster result into table\n";
				exit(1);
			}
		}
	}
	else {
		theMethod.fileinput(infile.c_str());
		theMethod.findCluster();
		theMethod.dumpTable("pcluster.tab", "cdsupdate.tab", cluster_id);
	}

	if (show) {
		theMethod.displayClusters();
	}
	cerr << "Total number of clusters: " << theMethod.getClusterCount() << endl;
	cerr << "Next cluster_id " << cluster_id << endl;

	return 0;
}

void usage(map<string, string> &param) {
	cerr << "Usage: orthopara\n" << "\t-h [host name default=" 
		<< param["host"] 
		<< "\n\t-d dbname default=" << param["dbname"]
		<< "\n\t-u user default=" << param["user"] 
		<< " for login to postgres server"
		<< "\n\t-i [cluster_id] starting number to make keys\n"
			  << "\t\tdefault=" << param["cluster_id"] << endl
		<< "\t-t [update_tab] the table containing the cds id\n"
			  << "\t\tdefault=" << param["update_tab"] << endl
		<< "\t-it [input_table] table for input\n"
			  << "\t\tdefault=" << param["input_tab"] << endl
		<< "\t-c [column name] column name to store cluster_id\n"
			  << "\t\tdefault=" << param["update_col"] << endl
		<< "\t-ct [cluster_tab] table name for storeing cluster result\n"
			  << "\t\tdefault=" << param["cluster_tab"] << endl
		<< "\t-cc [cluster_col] cluster column names\n"
			  << "\t\tdefault=" << param["cluster_col"] << endl
		<< "\t--show    show the clustering result to cout\n"
		<< "\t-q [query] query that produces the input for cluster\n"
			  << "\t\tdefault=" << param["query"] << endl
		<< "\t-f [input file] containing the relations\n";
	cerr << "\t-fo for file output\n";
	cerr << "\t-FO to turn off file output which is default\n";
}

