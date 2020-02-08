//#include "libpq++.h"
#include <string>
#include "dbinfo.h"
#include <fstream>
#include <iostream>
#include <pqxx>
#include <cstdlib>
#include <cstring>

using namespace std;
using namespace pqxx;

// file: picklongrep.cpp
// picks the longest representative from a group of same cds
/* the same_cds table is created as
 * create table same_cds(
 * 	cds_id integer primary key,
 * 	example integer not null);
 *
 * the query used for updating is 
 * select c.length, s.cds_id, s.example
 * from cds c join same_cds s on c.id=s.cds_id
 * order by s.example, c.length DESC, from_mRNA;
 * */

//void longrep(PgDatabase &db, string &of);
void longrep(result &qres, string &of);

int main(int argc, char* argv[]) {
	string intable="same_cds";
	if (argc < 2) {
		cerr << "Usage picklongrep -t same_cds\n";
		exit(1);
	}
	string ofile="same_cdslongrep.tab";
	bool FILE_OUT = true;
	int i=1;
	while (i<argc) {
		if (!strcmp(argv[i], "-t")) intable = argv[++i];
		else if (!strcmp(argv[i], "-of")) ofile = argv[++i];
		i++;
	}

	/* makes two connections to the database */
	////////////////////////////////////////////////////
	//PgDatabase pgdb("host=localhost dbname=orpara");
	connection pgdb("host=localhost dbname=orpara");
   /*
	if (pgdb.ConnectionBad()) {
		cerr << " connecting to pgdb failed\n";
		exit(1);
	}
   */
	//PgDatabase pgupdatedb("host=localhost dbname=orpara");
	connection pgupdatedb("host=localhost dbname=orpara");
   /*
	if (pgupdatedb.ConnectionBad()) {
		cerr << "Cannot open second db\n"; exit(1); 
	}
   */
	/////////////////////////////////////////////////////

	string query= "select i.example, c.id from ";
	query.append(intable).append(" i join cds c on c.id=i.cds_id order by i.example, c.length DESC, c.from_mrna");
	string updateQueryBase= "update ";
	updateQueryBase.append(intable).append(" set example='");

	string updateQuery;
	//const char *rep, *cdsid;
	string rep, cdsid;
	int updateCnt = 0;

	/* execute the ordering query for input */
   work Xaction(pgdb);
   result qres=Xaction.exec(query);
   /*
	if (!pgdb.ExecTuplesOk(query.c_str())) {
		cerr << "query\n" << query << "failed\n";
		exit(1);
	}
   */
	//cout << pgdb.Tuples() << " rows to update\n";
	cout << qres.size() << " rows to update\n";

	if (FILE_OUT) {
		longrep(qres, ofile);
		return 0;
	}

   work updateXaction(pgupdatedb);
			
	i=0;
	while (i<qres.size()) {
		if (i%200 == 0) cerr << ".";
		//rep = pgdb.GetValue(i, 0);     // example id in un-updated table
      qres[i][0].to(rep);
		//cdsid = pgdb.GetValue(i, 1);   // longest CDS id
      qres[i][1].to(cdsid);
		//if (strcmp(rep, cdsid)) {      // example is not the longest
		if (rep != cdsid) {      // example is not the longest
			updateCnt++;
			updateQuery=updateQueryBase;
			updateQuery.append(cdsid).append("' where example='");
			updateQuery.append(rep).append("'");

         updateXaction.exec(updateQuery);
         updateXaction.commit();
         /*
			if (!pgupdatedb.ExecCommandOk(updateQuery.c_str())) {
				cerr << updateQuery << " \nfailed!\n";
				exit(1);
			}
         */
		}
		i++;
      string tmp;
      qres[i][0].to(tmp);
		while (i<qres.size() && tmp == rep ) {
         i++;
         if (i<qres.size()) {
            qres[i][0].to(tmp);
         }
      }
	}
	cout << updateCnt << " rows updated\n";
	return 0;
}

/* db contains the query result */
void longrep(result &qres, string &of) {
	ofstream OU(of.c_str());

	int i=0;
	//const char *rep, *cdsid;
	string rep, cdsid;
	cerr << "total tuples=" << qres.size() << endl;
	while (i < qres.size()) {
		//rep = db.GetValue(i, 0);     
		qres[i][0].to(rep);
		//cdsid = db.GetValue(i, 1);   // longest CDS id
      qres[i][1].to(cdsid);
      string tmp;
      qres[i][0].to(tmp);
		while (i<qres.size() && tmp == rep) {
			OU << qres[i][1] << '\t' << cdsid << '\n';
			i++;
         if (i<qres.size()) qres[i][0].to(tmp);
		}
	}
	OU.close();
	cerr << "New long representative written to file: " << of << endl;
}
