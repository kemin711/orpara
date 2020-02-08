/** clusmethod.h
	 this is the database version of the orthopara clustering 
 	 protocol.  I used the integer as key.  String based key
	 will be developed later.  
	 It is simplified by taking PgDatabase as input for constructor.

	 Future developments should use file as input too
*/

#ifndef ORTHOMETHOD_H
#define ORTHOMETHOD_H
//#define HAVE_NAMESPACE_STD

//#include <multimap.h> in the std version, it is incuded in <map>
#include <map>
#include <vector>
//#include "libpq++.h"
#include <pqxx>
#include "cluster.h"

using namespace std;
using namespace pqxx;

/** remove all values with key k from mm and return the set
   */
set<int> grab(int k, multimap<int, int> &mm); 

class clusmethod
{
	public:
		/* connInfo is the connection info to the backend db, see
		 * lbpq definition.  Default is provided by a static
		 * char dbconnInfo[] in this class.
		 * The query must produce two columns of integers right now.
		 * seed->target both are keys
		 */
		//clusmethod(const char *query, PgDatabase &db);
		//clusmethod(const char *query, const char *connInfo = "host=localhost dbname=orpara");
		// replaced by simpler interface

		/** default constructor to make an empty object
		 * */
		clusmethod() { } 

		/**
		 * Version that does not depends on PgDatabase as a member, better
		 * Should be exactly like before.  The connection and query setting
		 * is done outside this class. 
		 * orthodb should contain the result of SQL with (query, target) as the
		 * tuple columns.
       *
       * not I have converted to the pqxx result type
		 */
		//clusmethod(PgDatabase &orthodb) { dbinput(orthodb); }
		clusmethod(result &orthodb) { dbinput(orthodb); }
		//void dbinput(PgDatabase &orthodb); 
      /** 
       * @param qres a query result
       * query result contain two columns
       * (obj1 ==> obj2)
       */
		void dbinput(result &qres); 
		/* take data from a flat file 
		 * first_number second_number
		 * use >> to read
		 * */
		clusmethod(char file[]) { fileinput(file); }
		void fileinput(const char file[]);

		/* this methods does all the work */
		vector<cluster>& findCluster();

		/* A function to display the two private members rel and revrel
		 * used for debuging.  */
		void printMap();

		/* Display the result to an output stream. Default to stdout */
		void displayClusters(ostream &os=cout);

		/* clears the input data to release memory */
		void clearMap() { rel.clear(); revrel.clear(); } 

		/* this methods loads the result into two database tables 
		 * update_tab = cds (update_col = pcluster)
		 * clutab = pcluster or ncluster (id, member) will be
		 * inserted
		 * id is the key for pcluster or ncluster
		 * cluster_id is the initial id number for naming clusters
		 * it will be incremented as each cluster is inserted.
		 * return true if successful
		 * */
		//bool loadTable(PgDatabase &orthodb, string update_tab, string update_col, 
		bool loadTable(connection &orthodb, string update_tab, string update_col, 
				string clutab, string clus_col, int &cluster_id);  
		void dumpTable(const char ins[], const char upd[], int &cluster_id); // write result into files

		int getClusterCount() { return allcl.size(); }

		/** remove all values with key k from mm and return the set
		 */
		//friend set<int> grab(int k, multimap<int, int> &mm); 

	private:
		// helper for findCluster
		set<int> processTargets(const set<int> &tset);  

		multimap<int, int> rel;    // forward relation
		multimap<int, int> revrel; // reverse relation
		//PgDatabase orthodb;  // just make things more complicated
		vector<cluster> allcl;  //store all the clusters; the result
};

#endif
