//#include "libpq++"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "group.h"
#include <cstdlib>

//int mapdiv(string div);
//bool trainpicker(const string &inf);  // return true if successful
bool goodQuality(gdiagnosis &ag, double zcut, ostream &sgl);
void usage(const string& infile, const string& configFile, double zcut, const string& m);
/* during debug stage, I output most of the data in human readable format
 * with diagnositic information.  At production stage, all data will be 
 * outputed as raw table format for data-loading into database.
 *
 * This program reads from a input files dumped from the maxdiv table.
 * */
int main(int argc, char* argv[])
{
	bool trainPicker = false;
	bool debug = true;
	string infile = "fr.train", configFile="/home/kzhou/proj/cluster/fr_pickdiv.istr";
	string matrixFile="/home/kzhou/proj/cluster/guid.mtx";
	double zcut = 1.2;
	int i = 1;
	while (i<argc) {
		if (!strcmp(argv[i], "-t")) trainPicker = true;
		else if (!strcmp(argv[i], "-c")) configFile = argv[++i];
		else if (!strcmp(argv[i], "-m")) matrixFile = argv[++i];
		else if (!strcmp(argv[i], "-z")) zcut = atof(argv[++i]);
		else if (!strcmp(argv[i], "-h")) { 
			usage(infile, configFile, zcut, matrixFile); return 1; 
		}
		else infile = argv[i];
		i++;
	}
	cerr << "zcut= " << zcut << endl;

	/* the gconst object must be first initialized to used the 
	 * child class, this action is essential for using both 
	 * gstat, and the group branch of objects
	 * */
	cerr << "gconst reading config file: " << configFile << " ...\n";
	gconst::readconf(configFile); 

	//////////// training  ////////////////////////
	if (trainPicker) {
		trainer theTrainer;
		if (theTrainer.train(infile)) return 0;
		else return 1;
	}
	ifstream IN(matrixFile.c_str());
	if (IN.fail()) {
		cerr << "openning file: " << matrixFile << " failed\n";
		exit(1);
	}
	gdiagnosis::readGuid(IN);
	IN.close();
	
	///// process true input data ///////////////
	//ifstream IN(infile.c_str());
	IN.open(infile.c_str());
	if (IN.fail()) { 
		cerr << "cannot open " << infile << "\n";
		exit(1);
	}
	char singlefile[] = "single.grp";  // verbose info
	char passfile[]="passed.grp"; 
	char passrawfile[] = "passed.raw";  // same as input
	//char failfile[] = "failed.ort";
	//char failrawfile[] = "failed.raw"; //raw file for computer processing
	char keyfile[] = "querydiv.tab"; // query, div passed only; passed+single
	// for database usage
	char editfile[] = "changed.grp";  // verbose info about how each is modified

	ofstream PASS(passfile);
	ofstream PASSRAW(passrawfile);
	//ofstream FAIL(failfile);
	//ofstream FAILRAW(failrawfile);
	ofstream SINGLE(singlefile);
	ofstream KEY(keyfile);   // for database loading
	ofstream ED(editfile);

	gdiagnosis ag(zcut);

	int query, passCnt=0, editCnt=0, singleCnt=0;
	bool moredata = true;
	IN >> query;
	cerr << "Processing data ...\n";
	while (moredata) {
		moredata = ag.next(query, IN);
tryagain :		if (ag.getDivCnt() == 1) { // group with single division
			ag.dump(SINGLE);
			ag.dumpKey(KEY);
			singleCnt++;
		}
		else if ( (ag.qcovpass() && ag.qualitypass()) ||
				(ag.qcovpass() && ag.scorepass()) ||
				(ag.scorepass() && ag.qualitypass())) { // best quality 1
			passCnt++;
			ag.dump(PASSRAW);  
			ag.dumpKey(KEY);
			if (debug) PASS << ag << "//////////////////////////\n\n";
		}
		else {
			gdiagnosis ag_copy(ag);  // for debuging
			vector<int> removedDiv = ag.rmlow(-3,3);
			if (removedDiv.empty()) { // nothing removed, group fine
				passCnt++;
				ag.dump(PASSRAW);
				ag.dumpKey(KEY);
				if (debug) PASS << ag << "/////////\n\n";
			}
			else {  // need to remove bad divisions
				if (debug) {  // for human to look at
					ED << ag_copy << endl << "Removed Divisions: ";
					for (int i=0; i<removedDiv.size(); i++) {
						ED << ag.divisions[removedDiv[i]] << " ";
					}
					ED << endl << "//////////\n\n";
				}
				editCnt++;
				goto tryagain;
			}
		}
	}
	IN.close();
	PASS.close();
	PASSRAW.close();
	//FAIL.close();
	KEY.close();
	ED.close();

	cerr << "Output written to passed: " << passfile 
		<< " single: " << singlefile << " key: " << keyfile << endl;

	cerr << "passed groups " << passCnt << "\tsingle groups " 
		<< singleCnt << "\tmodified groups " << editCnt << endl;
	cerr << "Use " << keyfile << " for loading database\n";

	return 0;
}

bool goodQuality(gdiagnosis &ag, double zcut, ostream &sgl) {
	if (ag.getDivCnt() > 1) {
		if (ag.isConserved() || ag.passedAvg(zcut) ) return true;
		else return false;
	}
	else ag.dump(sgl) << endl;
}

void usage(const string& infile, const string& configFile, double zcut, const string& m) {
	cerr << "Usage: pickortho infile default=" << infile << " if training\n"
		<< "\t-t do training of the program from the infile\n"
		<< "\t-c configFile default=" << configFile << endl
		<< "\t-z zcut default=" << zcut << endl
		<< "\t-m guidMatrix default=" << m << endl
		<< "\t-h print help message like this one\n";
}
