// cluster.h

#ifndef CLUSTER_H
#define CLUSTER_H

#include <set>
//#include <multimap.h>
#include <map>              // now map includes multimap
#include <string>
#include <iostream>
using namespace std;

class cluster 
{
	friend class clusmethod;
	public:
		cluster() {}
		cluster(int sd, const set<int> &tgs);

		/* constructs a new cluster from
		 * the first key->values, and removed
		 * the key and associated values from mm
		 */
		cluster(multimap<int, int> &mm); 
		void addSeed(int sd) { set1.insert(sd); }
		void addTargets(const set<int> &ts) { set2.insert(ts.begin(), ts.end()); }
		void display(ostream &os=cout);  // diagnostic function

		/* produce array in the format of postgres */
		string dumpArray() const;  
		/** convert to tab delimited text by default, or other delimiters
		 */
		string toText(const char delimiter='\t') const;

	 set<int> getElements() const;  //combine set1 and set2
	 int size() const { return set1.size() + set2.size(); }

	friend ostream& operator<<(ostream &os, const cluster &c);

	protected:
		set<int> set1;  // the seed set
		set<int> set2;  // the target set
};

#endif
inline cluster::cluster(int sd, const set<int> &tgs) {
	set1.insert(sd);
	set2 = tgs;
}

inline set<int> cluster::getElements() const { 
	set<int> tmp(set1);
	tmp.insert(set2.begin(), set2.end());
	return tmp;
}
	
