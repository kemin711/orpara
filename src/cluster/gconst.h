#ifndef GCONST_H
#define GCONST_H

#include <string>
#include <vector>

using namespace std;

/* this class is for constants that are used by group, and gstat class
 * this class is used as the super class, parent class
 * */

class gconst {
	public:
		////////// static public member for convinence /////////////
	static vector<string> divisions;

	static void readconf(const string &gdf);
	//static int getdividx(divstat* ds) { return getdividx(ds->id); }
	static int getdividx(const string &div) { return getdividx(div.c_str()); }
	static int getdividx(const char *div);
	// the following methods are never used, they are not usefull
	// at all
	static void setAnchorOrder(const int a[]);
	static void showAnchorOrder(ostream &ou);
	// returns the index of division based on index
	static int anchorIndex(const int dividx);

	///////////// static members //////////////////////////////
	static string pivotdiv;
	static int pivotdividx;     // the index of pivotdiv in divisions
	static string queryDiv;
	static float csvStd;         // conserved standard deviation 0.1
	static float csvStdHighVar;  // conserved std of highVarDiv 0.2
	static string highVarDiv[2]; 
	static float csvIdentity;    // if mean identity > csvIdentity, it is 
	static vector<int> searchOrder;  // anchor search order
};

#endif
