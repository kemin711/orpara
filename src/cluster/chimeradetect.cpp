#include "alnrange.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "libpq++.h"

using namespace std;

void checkcluster(const string &clid, PgDatabase *pgdb, int olpmg, const string &alntab);
void checkclusterFromStream(istream& ins);

// mark quality of aln result in alntab with -133 if 
// it is a weaker of the 2 or more chimeras
void checkchimera(PgDatabase *pgdb, int olpmg, const string &alntab, const string &matchtab);

// put result into an output stream ous
int ischimera(PgDatabase *pgdb, int olpmg, const string& matchtab);

enum prgtask { CHECKDB, QCDB };

int main(int argc, char* argv[]) {
	string clid;
	string infile;
	int olpmargin = 15;
	string alntab = "alntop";  // this is the current defaul
	// in the future alntop will be the default table
	if (argc<2) {
		cerr << "Usage chimeradetect -f inputfile of all half matches\n"
			<< " -c cluster_id\n" 
			<< " -m non_overlap margin default 15\n" 
			<< " --complete to detect all chimera from the halfaln table: alntophalfmatch\n"
			<< " -a alnsource_table default alntop\n"
			<< "--check inputtable\n";
		return 1;
	}

	string halfalntab;
	prgtask usr_request = CHECKDB;
	int i=1;
	while (i<argc) {
		if (!strcmp(argv[i], "-f")) infile = argv[++i];
		else if (!strcmp(argv[i], "-c")) clid = argv[++i];
		else if (!strcmp(argv[i], "-m")) olpmargin = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-a")) alntab = argv[++i];
		else if (!strcmp(argv[i], "--complete")) {
			halfalntab = argv[++i];
			usr_request = QCDB;
		}
		else if (!strcmp(argv[i], "--check")) {
			usr_request = CHECKDB;
			halfalntab = argv[++i];
		}
		i++;
	}
	// read input from file
	ifstream ins;
	if (!infile.empty()) {
		ins.open(infile.c_str());
		if (ins.fail()) {
			cerr << "Failed to open " << infile << endl;
			exit(1);
		}
		cerr << "File: " << infile << " opened for chimera checking!\n";
		checkclusterFromStream(ins);
		return 0;
	}

	PgDatabase pgdb("dbname=orpara host=localhost");
	if (!clid.empty())  {
		checkcluster(clid, &pgdb, olpmargin, alntab);
		return 0;
	}
	// fields: prt   len   mb me score identity ngidentity  cov
	// the order of the input is important: prt,mb,me is the correct one

	// will check every thing from the whole table
	if (!halfalntab.empty()) {
		if (usr_request == QCDB) {
			checkchimera(&pgdb,olpmargin, alntab, halfalntab);
		}
		else if (usr_request == CHECKDB) {
			cerr << "Checking all entries from table " 
				<< halfalntab << " for potetial chimeras\n";
			//string outfile="chimera.list";
			//ofstream outf(outfile.c_str());
			cerr << ischimera(&pgdb,olpmargin,halfalntab)
				<< " chimeras detected.\n";
		}
	}
	return 0;
}

/* Database base method
 * matchtab is the half match table
 * the output goes to alntab, it update the quality to -133 for this
 * sequence to all those that have lower scores if the
 * sequence is chimera
 * */
void checkchimera(PgDatabase *pgdb, int olpmg, const string &alntab, const string &matchtab) {
	//string matchtab = "alntophalfmatch";
	string query="select m.* from ( select prt from "
		+ matchtab + " where mb notnull group by prt having ( (max(mb)-min(mb)>70 and max(me)-min(me)>70 and max(mb)-min(me)>-40) or (stddev(mb)>70 and stddev(me)>70)) and count(*) > 2 except select prt from "
		+ matchtab + " where mb notnull and cov>0.9 and ngidentity>0.4) c join "
		+ matchtab + " m on c.prt=m.prt where m.mb notnull order by prt,mb,me";

	if (!pgdb->ExecTuplesOk(query.c_str())) {
		cerr << pgdb->ErrorMessage() << "\n" << query << "\nFailed\n";
		exit(1);
	}
	string seqid;
	int seqlen, i;
	//float score, identity, ngidentity, cov;

	int ntuples = pgdb->Tuples();
	if (ntuples == 0) {
		cerr << query << "\nProduce no result!\n";
		return;
	}
	cout << ntuples << " for candidate chimera\n";
	int r=0;
	while (r<ntuples) {
		seqid=pgdb->GetValue(r,0);
		seqlen=atoi(pgdb->GetValue(r,1));
		alnrange rr1(pgdb, r);

		vector<avgrange*> ars;
		ars.push_back(new avgrange(rr1));
		++r;

		while (r<ntuples && seqid == pgdb->GetValue(r,0)) {
			alnrange rr(pgdb, r);
			i=0;
			// the order of the testing is essential
			for (i=0; i<ars.size(); i++) {
				if (ars[i]->overlap(rr, olpmg)) {
					ars[i]->merge(rr);
					break;
				}
			}
			if (i == ars.size()) ars.push_back(new avgrange(rr));
			++r;
		}

		if (ars.size() == 1) {
			//cout << "--------------------------------\n";
			//cout << seqid << " len: " << seqlen << " is not chimera.\n";
			continue;
		}

		cout << "-------------------*******-------------\n";
		cout << "-- " << seqid << " " << seqlen << " is chimera.\n";

		// find the most significant match
		double maxscore=0;
		int maxidx;
		for (i=0; i<ars.size(); i++) {
			if (ars[i]->getScore() > maxscore) {
				maxscore=ars[i]->getScore();
				maxidx = i;
			}
		}
		for (i=0; i<ars.size(); i++) {
			if (i==maxidx) {
				cout << "\n-- Submatch retained.\n";
				ars[i]->sqlinfo(cout) << endl;
				continue;
			}
			cout << "\n-- Submatch removed.\n";
			ars[i]->sqlinfo(cout) << endl;
			// mark the quality of this as bad
			ostringstream ost;
			ost << "update " << alntab << " set quality=-133 where (query='"
				<< seqid << "' and qbegin >=" << ars[i]->minbegin()
				<< " and qend <= " << ars[i]->maxend() 
				<< ") or (target='" << seqid 
				<< "' and tbegin >= " << ars[i]->minbegin()
				<< " and tend <= " << ars[i]->maxend() << ")";
			string query=ost.str();
			cout << endl << query << ";\n\n";
			delete ars[i];
		}
	}
}

// the input table matchtab should be chimera_candidate, this
// can be changed to any table
// returns the number of chimera, -1 indicates failure
// put result into a default file chimera.list chimera.tab
int ischimera(PgDatabase *pgdb, int olpmg, const string& matchtab) 
{
	string query="select * from " + matchtab
			+ " order by prt, mb, me";

	if (!pgdb->ExecTuplesOk(query.c_str())) {
		cerr << pgdb->ErrorMessage() << "\n" << query << "\nFailed\n";
		exit(1);
	}
	string seqid;
	int seqlen, i;
	//float score, identity, ngidentity, cov;

	int ntuples = pgdb->Tuples();
	if (ntuples == 0) {
		cerr << query << "\nProduce no result!\n";
		return -1;
	}
	cout << ntuples << " for candidate chimera\n";
	string outfile="chimera.list";
	ofstream ous(outfile.c_str());
	ofstream tabous("chimera.tab");

	int r=0;
	int chimera_cnt=0;
	while (r<ntuples) {
		seqid=pgdb->GetValue(r,0);
		seqlen=atoi(pgdb->GetValue(r,1));
		alnrange rr1(pgdb, r);

		vector<avgrange*> ars;
		ars.push_back(new avgrange(rr1));
		++r;

		while (r<ntuples && seqid == pgdb->GetValue(r,0)) {
			alnrange rr(pgdb, r);
			i=0;
			// the order of the testing is essential
			for (i=0; i<ars.size(); i++) {
				if (ars[i]->overlap(rr, olpmg)) {
					ars[i]->merge(rr);
					break;
				}
			}
			// if not overlapping with existing range
			// add a new one
			if (i == ars.size()) ars.push_back(new avgrange(rr));
			++r;
		}

		// not a chimera
		if (ars.size() == 1) {
			//cout << "--------------------------------\n";
			//cout << seqid << " len: " << seqlen << " is not chimera.\n";
			continue;
		}

		chimera_cnt++;
		ous << "-------------------*******-------------\n";
		ous << "-- " << seqid << " " << seqlen << " is chimera.\n";

		// output result
		for (i=0; i<ars.size(); i++) {
			ous << *(ars[i]) << endl;
			tabous << seqid << "\t" << seqlen << "\t";
			ars[i]->writeTable(tabous);
			tabous << endl;
			delete ars[i];
		}
		ous << endl;
	}
	cerr << "Result written to chimera.list and chimera.tab\n";
	return chimera_cnt;
}

void checkclusterFromStream(istream& ins) {
	string prt, seqid;
	int len, mb, me, seqlen, i;
	float score, identity, ngidentity, cov;

	ins >> prt >> len >> mb >> me >> score >> identity >> ngidentity >> cov;
	while (!ins.eof()) {
		alnrange rr1(mb,me,score,ngidentity,cov);
		vector<avgrange*> ars;
		ars.push_back(new avgrange(rr1));
		seqid = prt;
		seqlen = len;
		ins >> prt >> len >> mb >> me >> score >> identity >> ngidentity >> cov;

		// working on the same sequence
		while (!ins.eof() && seqid == prt) {
			alnrange rr(mb,me,score,ngidentity,cov);

			i=0;
			// the order of the testing is essential
			for (i=0; i<ars.size(); i++) {
				if (ars[i]->overlap(rr)) {
					ars[i]->merge(rr);
					break;
				}
			}
			if (i == ars.size()) ars.push_back(new avgrange(rr));
			ins >> prt >> len >> mb >> me >> score >> identity >> ngidentity >> cov;
		}
		if (ars.size() == 1) {
			cout << seqid << " " << seqlen << " is not chimera.\n";
		}
		else {
			cout << seqid << " " << seqlen << " is chimera.\n";
		}

		for (i=0; i<ars.size(); i++) {
			cout << *ars[i] << endl;
			delete ars[i];
		}
	}
}

// based on postgres database, work on one cluster
void checkcluster(const string &clid, PgDatabase *pgdb, int olpmg, const string &alntab) {
	string matchtab = "clmatch";
	if (alntab == "alntop") matchtab = "seqmatch";

	string query="select m.* from (select prt from "
		+ matchtab + "(" + clid 
		+ ") group by prt having (max(mb)-min(mb)>70 and max(me)-min(me)>70 and max(mb)-min(me)>-30) or (stddev(mb)>70 and stddev(me)>70) except select prt from "
		+ matchtab + "(" + clid + ") where cov>0.9 and ngidentity>0.4) c join "
		+ matchtab + "(" + clid + ") m on c.prt=m.prt";

	if (!pgdb->ExecTuplesOk(query.c_str())) {
		cerr << pgdb->ErrorMessage() << "\n" << query << "\nFailed\n";
		exit(1);
	}
	string seqid;
	int seqlen, i;
	//float score, identity, ngidentity, cov;

	int ntuples = pgdb->Tuples();
	if (ntuples == 0) {
		cerr << query << "\nProduce no result!\n";
		return;
	}
	cout << ntuples << " for candidate chimera\n";
	int r=0;
	while (r<ntuples) {
		seqid=pgdb->GetValue(r,0);
		seqlen=atoi(pgdb->GetValue(r,1));
		alnrange rr1(pgdb, r);

		vector<avgrange*> ars;
		ars.push_back(new avgrange(rr1));
		++r;

		while (r<ntuples && seqid == pgdb->GetValue(r,0)) {
			alnrange rr(pgdb, r);
			i=0;
			// the order of the testing is essential
			for (i=0; i<ars.size(); i++) {
				if (ars[i]->overlap(rr, olpmg)) {
					ars[i]->merge(rr);
					break;
				}
			}
			if (i == ars.size()) ars.push_back(new avgrange(rr));
			++r;
		}

		if (ars.size() == 1) {
			//cout << "--------------------------------\n";
			//cout << seqid << " len: " << seqlen << " is not chimera.\n";
			continue;
		}

		cout << "-------------------*******-------------\n";
		cout << "-- " << seqid << " " << seqlen << " is chimera.\n";

		// find the most significant match
		double maxscore=0;
		int maxidx;
		for (i=0; i<ars.size(); i++) {
			if (ars[i]->getScore() > maxscore) {
				maxscore=ars[i]->getScore();
				maxidx = i;
			}
		}
		for (i=0; i<ars.size(); i++) {
			if (i==maxidx) {
				cout << "\n-- Submatch retained.\n" << *ars[i] << endl;
				continue;
			}
			cout << "\n-- Submatch removed.\n"  << *ars[i] << endl;
			// mark the quality of this as bad
			ostringstream ost;
			ost << "update " << alntab << " set quality=-133 where (query='"
				<< seqid << "' and qbegin >=" << ars[i]->minbegin()
				<< " and qend <= " << ars[i]->maxend() 
				<< ") or (target='" << seqid 
				<< "' and tbegin >= " << ars[i]->minbegin()
				<< " and tend <= " << ars[i]->maxend() << ")";
			string query=ost.str();
			cout << endl << query << ";\n\n";
			delete ars[i];
		}
	}
}

