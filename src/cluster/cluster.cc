// cluster.cc
// using the id, integer as key

#include "cluster.h"
#include <sstream>  // do not use strstream in the old C++
#include <iterator>

cluster::cluster(multimap<int, int> &mm) {
	multimap<int, int>::iterator p = mm.begin();
	int firstKey = p->first;

	set1.insert(firstKey);
	while (p != mm.upper_bound(firstKey)) set2.insert(p++->second);
	mm.erase(firstKey);
}

void cluster::display(ostream &os) {
	set<int>::iterator p = set1.begin();
	os << "seed: ";
	while (p != set1.end()) {
		os << *p << " ";
		p++;
	}
	os << "==> target: ";
	p = set2.begin();
	while (p != set2.end()) {
		os << *p << " ";
		p++;
	}
	os << endl;
}

ostream& operator<<(ostream &os, const cluster &c) {
	copy(c.set1.begin(), c.set1.end(), ostream_iterator<int>(os, " "));
	os << "==> ";
	copy(c.set2.begin(), c.set2.end(), ostream_iterator<int>(os, " "));
	return os;
}

string cluster::toText(const char delimiter) const {
	ostringstream ostr;  // this object got destroyed when out of scope
	set<int>::const_iterator si;
	for (si=set1.begin(); si != set1.end(); si++) ostr << *si << delimiter;
	si=set2.begin();
	ostr << *si++;
	while (si != set2.end()) {
		ostr << delimiter << *si;
		si++;
	}
	return ostr.str();
}


string cluster::dumpArray() const {
	ostringstream ostr;  // this object got destroyed when out of scope
	ostr << '{';
	set<int>::const_iterator si;
	for (si=set1.begin(); si != set1.end(); si++) ostr << *si << ',';
	si=set2.begin();
	ostr << *si++;
	while (si != set2.end()) {
		ostr << ',' << *si;
		si++;
	}
	ostr << '}';
	return ostr.str();
}

