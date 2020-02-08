#ifndef SCOREPAIR_H
#define SCOREPAIR_H
//#define HAVE_NAMESPACE_STD

#include <string>
#include <map>
#include <vector>
#include <set>
//#include "forest.h"
#include "hatrees.h"
//#include "libpq++.h"
#include <pqxx>

using namespace std;
using namespace pqxx;

//enum pairquality { good, bad, modify };
//void deleteCluster(PgDatabase &pgdb, const string &clusterid);
void deleteCluster(connection &pgdb, const string &clusterid);

class scorepair {
	public:
		void read(istream &in, int cut);
		void clear() { names.clear(); seqlens.clear(); rel.clear(); }

		// 1=good, 9=bad, 5=modify; reading subclstr
		int quality();
		void showcluster(ostream &ou);

		/* return all those that are not included in any cluster */
		set<string> badMembers() const;
		//void loadTable(PgDatabase &pgdb);
		void clusterArray(vector<set<string> > &clar);
		vector<set<string> > getClusters();  // too expensive don't call

		void addNewCluster(connection &pgdb, const string &clusterid, 
				int &clusterIdCnt);

	private:
		multimap<string, string> rel;  // use integer for speed
		vector<string> names;  // name -> identifier
		vector<int> seqlens;
		hatrees<string> subclstr;
};

#endif
