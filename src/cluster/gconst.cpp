#include "gconst.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

vector<string> gconst::divisions=vector<string>();
string gconst::pivotdiv = "human";
int gconst::pivotdividx = 0;

vector<int> gconst::searchOrder = vector<int>();
string gconst::queryDiv = "fugu";  // will be overwritten
string gconst::highVarDiv[2] = {"bony", "vrt"};
float gconst::csvStd=0.1;
float gconst::csvStdHighVar=0.2;
float gconst::csvIdentity=0.85;

void gconst::setAnchorOrder(const int a[]) {
	for (int i=0; i<searchOrder.size(); i++) searchOrder[i] = a[i];
}
void gconst::showAnchorOrder(ostream &ou) {
	for (int i=0; i<searchOrder.size(); i++) ou << searchOrder[i] << '\t';
}

int gconst::anchorIndex(const int dividx) {
	for (int i=0; i<searchOrder.size(); i++) {
		if (dividx == searchOrder[i]) return i;
	}
	return -1;  // anchor not found, out of ranges, > divisions.size()
}
int gconst::getdividx(const char* div) {
	int i = 0;
	while (i<divisions.size()) {
		if (div == divisions[i]) return i;
		i++;
	}
	cerr << div << " not found in our divisions\n";
	exit(1);
}
													    
/* reading guid file, instruction file */
//void group::readconf(const char gdf[]) {
void gconst::readconf(const string &gdf) {
	ifstream GUID(gdf.c_str());
	if (GUID.fail()) { cerr << " guid file open failed\n"; exit(1); }
	string tmp;
	GUID >> tmp;
	while (tmp[0]=='/' and tmp[1]=='/') {
		getline(GUID, tmp);  // discarding the rest of the line
		GUID >> tmp;
	}
	if (tmp != "QUERY") { 
		cerr << "The first non-comment line should be QUERY\n";
		exit(1);
	}
	GUID >> queryDiv;
	GUID >> tmp;
	while (tmp[0]=='/' and tmp[1]=='/') {
		getline(GUID, tmp);
		GUID >> tmp;
	}
	if (tmp != "TARGETS") {
		cerr << "There should be a target line\n";
		exit(1);
	}
	GUID >> tmp;
	while (tmp != "SEARCHORDER") {
		divisions.push_back(tmp);
		GUID >> tmp;
	}
	int xx;
	for (int i=0; i<divisions.size(); i++) {
		GUID >> xx;
		searchOrder.push_back(xx);
	}
	GUID >> tmp;
	while (tmp != "CONSERVED") {
		getline(GUID, tmp);
		GUID >> tmp;
	}
	GUID >> csvIdentity >> csvStd >> csvStdHighVar >> highVarDiv[0] >> highVarDiv[1];
	GUID >> tmp;
	while (tmp != "PIVOTDIV") {
		getline(GUID, tmp);
		GUID >> tmp;
	}
	GUID >> pivotdiv;
	GUID >> tmp;
	while (tmp != "PIVOTDIV_IDX") {
		getline(GUID, tmp);
		GUID >> tmp;
	}
	GUID >> pivotdividx;
}

