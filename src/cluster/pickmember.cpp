//#include "libpq++"
#include <iostream>
#include <fstream>
#include <string>

int mapdiv(string div);

int main(int argc, char* argv[])
{
	ifstream IN("rel.tab");
	if (IN.fail()) { 
		cerr << "cannot open rel.tab\n";
		exit(1);
	}
	int query, target, matchlen;
	string div[4];
	float identity, qcov, tcov, score;
	double avgscore[8], onescore[8];

	IN >> query >> div >> target >> target >> identity >> matchlen >> qcov
		>> tcov >> score;
	while (!IN.eof()) {
		int currentQuery = query;
		while (currentQuery == query) {
			onescore[mapdiv(div)] = score;

	return 0;
}

int mapdiv(string div) {
	if (div == "amp") return 0;
	else if (div == "fur") return 1;
	else if (div == "mam") return 2;
	else if (div == "mum") return 3;
	else if (div == "pri") return 4;
	else if (div == "rod") return 5;
	else if (div == "sau") return 6;
	else if (div == "vrt") return 7;
}
