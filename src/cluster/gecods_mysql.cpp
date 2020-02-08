// file: gecods_mysql.cpp
#include "hatrees.h"
//#include <stack>
#include <map>
//#include "/usr/local/pgsql/include/libpq++.h"
//#include "libpq++.h"
#include <mysql++.h>
#include <string>
//#include "loadpgdb.h"
//#include <multimap.h>
#include <iostream>

using namespace std;
using namespace mysqlpp;

/** 
 * This program is adapted to work with the mysql C++ driver.
 * It is based on the hatrees library.
 * */

void usage(ostream &os);
enum IdType {STRID, INTID};

void createOutputTable(Connection &conn, const string &outab, IdType idtype);
template <class T>
void loadIntoTable(hatrees<T> &cl, const string &outab, Connection &conn, IdType idtype);
template<class T> 
void clusterFromFile(hatrees<T> &cluster, ofstream &OUT, const string &infile);
template<class T>
void clusterFromTable(UseQueryResult &res, hatrees<T> &clusters);

int main(int argc, char* argv[]) {
   string host, database; 
   string user="kzhou";
   string password="Hotcreek2800";

	string infile="clinput.tab"; 
	string input_table="link"; // default input_table name
	string output_table="cluster";
	string output_file="cluster.tab";
	bool fileInput = false;
	bool fileOutput = true;  // a lot faster for loading
	string inputQuery;
   /* default idtype set to string which is works in
    * every case */
	IdType idtype = STRID; // 0 for string as id, 1 for integer as id
	int clid = 0; // the starting number for clusterid which is auto_increment
	
	if (argc < 2) usage(cerr);
	int i=1;
	while (i<argc) {
		if (!strcmp(argv[i], "-f")) {
			infile=argv[++i];
			fileInput = true;
			fileOutput = true;
		}
		else if (!strcmp(argv[i], "-h")) host=argv[++i];
		else if (!strcmp(argv[i], "-d")) database=argv[++i];
		else if (!strcmp(argv[i], "-i")) input_table=argv[++i];
		else if (!strcmp(argv[i], "--idtype") || !strcmp(argv[i], "-k")) {
			++i;
			if (!strcmp(argv[i],"int") || !strcmp(argv[i], "INT")) {
				idtype = INTID;
			}
			else if (!strcmp(argv[i], "str") || !strcmp(argv[i], "STR")
               || !strcmp(argv[i], "string")) 
            idtype = STRID;
			else {
				cerr << "idtype must be either int or str or string\n";
				exit(1);
			}
		}
		else if (!strcmp(argv[i], "--clid-serial")) {
			clid = atoi(argv[++i]); // initial cluster id
         if (clid<1) { 
            cerr << "auto_increment id must start with 1 or greater\n";
            exit(1);
         }
		}
		else if (!strcmp(argv[i], "--clid-rep")) {
         clid = 0;
      }
		else if (!strcmp(argv[i], "-o")) { 
			fileOutput = false;
			output_table=argv[++i];
		}
		else if (!strcmp(argv[i], "-of")) {
			if (++i >= argc) {
				cerr << "you must supply a fileName after -of option\n";
				exit(1);
			}
			output_file=argv[i];
		}
		else if (!strcmp(argv[i], "-q")) inputQuery = argv[++i];
		else { 
			cerr << argv[i] << " is not a legal option\n";
			usage(cerr); 
		}
		i++;
	}

   // for file output
	hatrees<string> clusters;
	hatrees<unsigned int> int_clusters;
	ofstream OUT;
	if (fileOutput) {
		OUT.open(output_file.c_str());
		if (OUT.fail()) { exit(1); }
	}
   // for file as input
	if (fileInput) {
      if (idtype == STRID) clusterFromFile(clusters, OUT, infile);
      else clusterFromFile(int_clusters, OUT, infile);
		return 0;
	}

	////// working from a database /////////////////////
   //onst char *db, const char *host="", const char *user="", const char
   //*passwd="", uint port=0, my_bool compress=0, unsigned int
   //connect_timeout=60, cchar *socket_name=0, unsigned int client_flag=0)
   try {
      cerr << "Trying to connect to " << database << endl;
      Connection conn(database.c_str(), host.c_str(), user.c_str(), password.c_str());
      if (conn.connected()) {
         cerr << "Connected to " << database << " on "
            << host << endl;
      }
      else {
         cerr << "Failed to connect to " << database << endl;
         return 1;
      }
      Query query=conn.query();
      if (inputQuery.empty()) inputQuery = "select * from " + input_table;
      query << inputQuery;
      UseQueryResult res = query.use();
      if (idtype == STRID) {
         clusterFromTable(res, clusters);
         if (fileOutput) {
            if (clid == 0) clusters.showCluster(OUT); 
            else clusters.showClusterIntId(OUT, clid);
         }
         else loadIntoTable(clusters, output_table, conn, idtype);
      }
      else { // int id type
         clusterFromTable(res, int_clusters);
         if (fileOutput) {
            if (clid == 0) int_clusters.showCluster(OUT);
            else int_clusters.showClusterIntId(OUT, clid);
         }
         else loadIntoTable(int_clusters, output_table, conn, idtype);
      }
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      return -1;
   }

	return 0;
}

void usage(ostream &os) {
	os << "Usage: gecods -f infile\n"
		<< "\t-d db           required\n"
		<< "\t-h host         required\n"
      << "\t-f infile       required\n"
		<< "\t-i input_table  required infile or input_table\n"
		<< "\t-o output_table default=cluster\n"
		<< "\t-of output_file will only output to file\n"
		<< "\t-q SQL_query    default=select query,target from input_table\n"
		<< "\t-k key_type     int or string  default=string\n"
		<< "\t--idtype [str, int] if --idtype set to int then needs clid\n"
		<< "\t--clid-serial integer >=0 initial cluster id, default=1\n"
		<< "\t--clid-rep     flag to pick one member's id at random as the cluster id\n"
		<< "\t\toptional\n";
	os << "file input or db input are exclusive\nIf -i is used then the program will take input from a table. The first two columns must contains the input data.\n";
	exit(1);
}

template<class T>
void clusterFromTable(UseQueryResult &res, hatrees<T> &clusters) {
   multimap<T,T> pairs;
   res.disable_exceptions();
   Row row = res.fetch_row();
   while (row) {
      pairs.insert(make_pair(static_cast<T>(row.at(0)), static_cast<T>(row[1])));
      row=res.fetch_row();
   }
   cerr << pairs.size() << " Input stored, started to do the cluster\n";
   clusters.readFromMap(pairs);
   cout << clusters.getNodeCount() << " unique members\n";
}

template<class T> 
void clusterFromFile(hatrees<T> &cluster, ofstream &OUT, const string &infile) {
   cluster.readFromFile(infile);
   OUT << "RESULTS: " << cluster.getNodeCount() << "\n";
   cluster.showStore(OUT);
   cout << "Now it is sorted according to cluster id\n";
   cluster.showCluster(OUT);
   cluster.showClusterByLine(cout);
}

void createOutputTable(Connection &conn, const string &outab, IdType idtype) {
   Query query=conn.query();
   try {
      query.exec("drop table if exists " + outab);
      if (idtype == STRID) 
         query.exec("create table " + outab + " (rep varchar(50), member varchar(50) primary key, key r(rep))");
      else 
         query.exec("create table " + outab + " (rep integer unsigned, member integer unsigned primary key, key r(rep))");
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
}

/* insert into a table, use the rep's id as the cluster id.
 * I have not written one version with auto_increment integer id yet.
 * */
template<class T>
void loadIntoTable(hatrees<T> &cl, const string &outab, Connection &conn, IdType idtype) {
   int commitInterval = 1000;
   bool firstValue=true;
   createOutputTable(conn, outab, idtype);
   try {
      Query query=conn.query();
      query << "insert delayed into " << outab << " values ";
      int i=0;
      multimap<T,T> clresult = cl.getCluster();
      typename multimap<T,T>::const_iterator it=clresult.begin();

      while (it != clresult.end()) {
         if (firstValue) {
            firstValue=false;
         }
         else {
            query << ",";
         }
         query << "('" << it->first << "', '" << it->second << "')";
         ++i;
         if (i % commitInterval == 0) {
            query.execute();
            query.reset();
            query << "insert delayed into " << outab << " values ";
            firstValue=true;
         }
         ++it;
      }
      if (i % commitInterval != 0) {
         query.execute();
      }
      cerr << clresult.size() << " Result stored in table: " << outab << endl;
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
}
