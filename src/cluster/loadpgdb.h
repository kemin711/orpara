#ifndef LOADPGDB_H
#define LOADPGDB_H
/* commonly used routines for loading results into database */
//#define HAVE_CXX_STRING_HEADER

#include <string>
#include <map>
#include <vector>
#include <set>
//#include "forest.h"
//#include "libpq++.h"
#include <pqxx>


void deleteCluster(pqxx::connection &pgdb, const std::string &clusterid);
//void addNewCluster(PgDatabase &pgdb, const string &clusterid, int &clusterIdCnt);
void loadClusterIntoDB(pqxx::connection &pgdb, 
		const std::multimap<int, int> &clusters, int &clusterIdCnt);

#endif
