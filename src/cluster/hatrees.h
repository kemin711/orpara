// file hatrees.h  replaces forest.h
// template class;  no *.cpp file associated with it
#ifndef DISET_H
#define DISET_H
//#define HAVE_NAMESPACE_STD

#include <string>
#include <map>
//#include <hash_map.h>
//#include <multimap.h> included by <map>
#include <utility>   // contains pair.h
#include <fstream>
#include <sstream>
//#include "libpq++.h"
//#include <pqxx>
#include <set>
#include <vector>
#include <iostream>
#include <cstdlib>

//#include <istream.h>
//
/** hash_map behave abnormally under template class and functions
 * but map (tree_map) works for this algorithm 
 *
 * I am getting rid of the libpqxx dependency.
 * This requires that postgress being installed on the 
 * system. This is only an convinience not a necessity.
 * This is a header only class.
 * */

//#define USE_HASH_MAP  // will insert duplicate key

// there are a lot std utilites used in this 
// header, the following statement polute the namespace
using namespace std;
//using namespace pqxx;

template<class T> struct node
{
	node() :  key(), parent(this), rank(0)  { }  // point to self
	node(const T &k) : key(k), parent(this), rank(0) { }
	node(const T &k, node* p) : key(k), parent(p), rank(0) { }
	/* set parent = this is essential */
	node(const node &n) : key(n.key), parent(this), rank(n.rank) { }
	// if use parent = n.parent then segmentation fault
	~node() { /*nothing needs to be done*/ }
	bool operator==(const node& n) const { return key==n.key; }
	bool operator<(const node &n) const { return key < n.key; }
	//node& operator=(const node &n) { if (&n != this) { key = n.key; parent = n.parent; rank = n.rank; } return *this; }
	node& operator=(const node &n) { if (&n != this) { key = n.key; parent = this; rank = n.rank; } return *this; }

	T key;   // stored in the map<key, value>
	node<T>* parent;
	int rank;
};

/* tree operations */
///////////////////////////////////////////////////////
/* x is not altered */
template<class T> node<T>* getRoot(node<T> *x) {
	if (x != x->parent)  x->parent = getRoot(x->parent); 
	return x->parent;
}
template<class T> void mergeRoot(node<T> *&r1, node<T> *&r2) {
	if(r1 == r2) return;  // if both are the same set
	if (r1->rank > r2->rank) r2->parent = r1;
	else {
		r1->parent = r2;
		if (r1->rank == r2->rank) r2->rank++;
	}
}
template<class T> void join(node<T> &x, node<T> &y) {
	node<T> *r1 = getRoot(&x);  // this is necessary for the compiler
	node<T> *r2 = getRoot(&y);
	mergeRoot(r1, r2);
	//mergeRoot(getRoot(&x), getRoot(&y));  // compiler got confused 
}
template<class T> void join(node<T>* x, node<T>* y) {
	node<T> *r1 = getRoot(x);  // this is necessary for the compiler
	node<T> *r2 = getRoot(y);
	mergeRoot(r1, r2);
	//mergeRoot(getRoot(x), getRoot(y));  // compiler got confused 
}

////////////////// hash_map specification /////////////////////////////
/* for hash_map of string as key, hash is not working
 * should be implemented later */
/*
template<> struct hash<string> {
	size_t operator()(const string &key) const {
		size_t res=0;
		typedef string::const_iterator CI;
		CI  p = key.begin();
		CI end= key.end();

		while (*p) {
			res = (res<<4) + *p++;
			size_t g = res & 0xF0000000L;
			if (g) res ^= g >> 24;
			res &= ~g;
		}
		return res;
	}
};
*/
///////////////////////////////////////////////////////
//   hatrees class definition /////////////////////////
/** 
 * this class provides simple file output methods.
 * No database table output is defined in this class
 * bedause it will require different database drivers.
 * The getCluster() method will produce a relational table.
 */

template<class T> class hatrees
{
	public:
		hatrees() { }
      /** build a cluster object directly from multimap
       * Save the call to readFromMap()
       */
		hatrees(const std::multimap<T,T> &mm) { readFromMap(mm); }
		~hatrees();
		void clear() { nodes.clear(); result.clear(); }
		///////// input operations /////////////////
		void readFromMap(const std::multimap<T, T> &m);
		void readFromFile(const string &file);
      /* only gecods.cpp needs this function
       * I am not compiling gecods.cpp anymore
       */
		//void readFromDB(PgDatabase &db);
		//void readFromDB(pqxx::result &qres);

		/********** output ********************/
		/** using representative as cluster id, results inserted into tab */
		//void loadTable(PgDatabase &db, string tab);
		//void loadTable(pqxx::connection &db, string tab);

		/** changes the parent pointer of the nodes	*/
		void showStore(ostream &os) const;  

		/** number of unique members */
		int getNodeCount() const { return nodes.size(); }
		/** 
       * in a table format: (cluster_id(representative), members) */
		void showCluster(ostream &os, bool reverse= true);
		/** an initial id is given, the program will
       * auto_increment this id. This should work with
       * either string or int type of idtype.
       * */
		void showClusterIntId(ostream &ous, int &id);

		/** return the number of clusters */
		int showClusterByLine(ostream &os);
      /** return all the keys as a set<T> it is all the 
       * members in the input. Useful for finding singleton
       * clusters. */
		set<T> keyset() const;
		vector<T> keyarray() const;

		/** produce an array of clusters, used call by reference
		 * to avoid creating tmp and copying.
       * The result could be used for further analysis. 
       * This is one of the output methods.
       * */
		void clusterArray(vector<set<T> > &vecset);

      /** Return a relational table
       * with two columns: | representative | member |
       * Use this method to load into database tables.
       * Singletons are not included in the output.
       */
		std::multimap<T,T>& getCluster() { if (result.empty()) transform(); return result; }

#ifdef USE_HASH_MAP
	//typename hash_map<T, node<T>* >::iterator niterator;
	//typename hash_map<T, node<T>* >::const_iterator const_niterator;
	typedef typename hash_map<T, node<T>* >::iterator niterator;
	typedef typename hash_map<T, node<T>* >::const_iterator const_niterator;
#else
	typedef typename map<T, node<T>* >::iterator niterator;
	typedef typename map<T, node<T>* >::const_iterator const_niterator;
	//typename map<T, node<T>* >::iterator niterator;
	//typename map<T, node<T>* >::const_iterator const_niterator;
#endif

	private:
#ifdef USE_HASH_MAP 
		hash_map<T, node<T>* > nodes;
#else
		map<T, node<T>* > nodes;
#endif
		/** This is the result table:
       * representative->members 
       * Should have another method to return a reversed result.*/
		//multimap<T, T> cluster;  // sorted result, loaded after
		std::multimap<T, T> result;  // sorted result in map format

		/* transform the nodStore cluster into mutimaped cluster
		 * in a table format cluster_id -> members */
		void transform();
};

////////// hatrees class template methods //////////////////////////
template<class T> hatrees<T>::~hatrees() { 
	niterator MI = nodes.begin();
	while (MI != nodes.end()) {
		delete MI->second;
		MI++;
	}
}

template<class T> void hatrees<T>::readFromMap(const std::multimap<T, T> &m) {
	pair<niterator, bool> p1, p2;
	T f1, f2; 
	typename std::multimap<T, T>::const_iterator imi = m.begin();

	while (imi != m.end()) {
		f1 = imi->first;
		f2 = imi->second;
		p1 = nodes.insert(pair<T, node<T>* >(f1, new node<T>(f1)));
		p2 = nodes.insert(pair<T, node<T>* >(f2, new node<T>(f2)));
		join(p1.first->second, p2.first->second);
		imi++;
	}
}

template<class T> void hatrees<T>::readFromFile(const string &file) {
	ifstream IN(file.c_str());
	if (IN.fail()) {
		cerr << "reading from " << file << "failed\n";
		exit(1);
	}
	string ln;
	pair<niterator, bool> p1, p2;
	T f1, f2; 
	//hash_map<T, T> xx;   // hash_map inserts duplicate key under template
	//hash_map<T, int> hm;
	//map<T, T> yy;
	//int id=0;

	getline(IN, ln);
	while (!IN.eof()) {
		istringstream ist(ln);
		ist >> f1 >> f2;
		/*
		xx.insert(pair<T, T>(f1, f2));
		yy.insert(pair<T, T>(f1, f2));
		if (hm.find(f1) == hm.end()) {
			hm.insert(pair<T, int>(f1, id++));
		}
		if (hm.find(f2) == hm.end()) {
			hm.insert(pair<T, int>(f2, id++));
		}
		*/
		p1 = nodes.insert(pair<T, node<T>* >(f1, new node<T>(f1)));
		p2 = nodes.insert(pair<T, node<T>* >(f2, new node<T>(f2)));
		join(p1.first->second, p2.first->second);
		getline(IN, ln);
	}
	/*
	cout << "Size of hash_map<string, string> is " << xx.size() << endl;
	cout << "Size of map<string, string> is " << yy.size() << endl;
	cout << "Size of hash_map<string, int> is " << hm.size() << endl;
	hash_map<T, T>::iterator xxi = xx.begin();
	map<T,T>::iterator yyi = yy.begin();
	while (xxi != xx.end()) {
		cout << xxi->first << "\t" << xxi->second << endl;
		xxi++;
	}
	cout << "================================\n";
	while (yyi != yy.end()) {
		cout << yyi->first << "\t" << yyi->second << endl;
		yyi++;
	}
	cout << "-----------------------------\n";
	hash_map<T, int>::iterator hmi = hm.begin();
	while (hmi != hm.end()) {
		cout << hmi->first << "\t" << hmi->second << endl;
		hmi++;
	}
	*/
}

/* db has the result ready; will use field 1, and 2 as input */
//template<class T> void hatrees<T>::readFromDB(PgDatabase &db) {
/* removing this function so that this library will not
 * dependent on postgres
template<class T> void hatrees<T>::readFromDB(pqxx::result &qres) {
	pair<niterator, bool> p1, p2;
	T f1, f2; 
	string ln, tmp;

	int i;
	for (i=0; i< qres.size(); i++) {
		//ln = db.GetValue(i,0);
		qres[i][0].to(ln);
      qres[i][1].to(tmp);
		//ln.append("\t").append(db.GetValue(i,1));
		ln.append("\t").append(tmp);
      
		istringstream ist(ln);
		ist >> f1 >> f2;
		p1 = nodes.insert(pair<T, node<T>* >(f1, new node<T>(f1)));
		p2 = nodes.insert(pair<T, node<T>* >(f2, new node<T>(f2)));
		join(p1.first->second, p2.first->second);
	}
}
*/
/* load into relational table(member_key, member_representative) */
//template<class T> void hatrees<T>::loadTable(PgDatabase &db,string tab) {
/* disable to make this library less dependent on postgress 
 * installation
template<class T> void hatrees<T>::loadTable(pqxx::connection &db,string tab) {
	string queryBase = "insert into " + tab + " values(";
	
	niterator it = nodes.begin();
	node<T>* root;
	cerr << "Loading into table member_key\tgroup_representative\n";
   work Xaction(db, "forinsertion");

	while (it != nodes.end()) {
		root = getRoot(it->second->parent);
		string query = queryBase;
		ostringstream ous;
		ous << '\'' << it->first << "','" << root->key << "')";
		query += ous.str();
      try {
         Xaction.exec(query);
      }
      catch (exception &err) {
         cerr << err.what() << endl;
         exit(1);
      }
      Xaction.commit();
      */
      /*
		if (!db.ExecCommandOk(query.c_str())) {
			cerr << query << " failed\n";
			exit(1);
		}
      */
/*
		it++;
	}
	cerr << "Total number of unique keys is " << nodes.size() << endl;
}
*/

/* members->cluster_id sorted according to members if using <map>
 * if using hash_map then not in any order */
template<class T> void hatrees<T>::showStore(ostream &os) const {
	const_niterator it = nodes.begin();
	node<T>* root;
	os << "key\tcluster_representative\n";
	while (it != nodes.end()) {
		root = getRoot(it->second->parent);
		os << it->first << "\t" << root->key << endl;
		it++;
	}
	os << "Total number of unique keys is " << nodes.size() << endl;
	cerr << "Total number of unique keys is " << nodes.size() << endl;
}

/* return the key as a set<T> */
template<class T> set<T> hatrees<T>::keyset() const {
	set<T> tmpset;
	const_niterator it = nodes.begin();
	while (it != nodes.end()) {
		tmpset.insert(it->first);
		it++;
	}
	return tmpset;
}

template<class T> vector<T> hatrees<T>::keyarray()  const {
	vector<T> tmp;
	const_niterator it = nodes.begin();
	while (it != nodes.end()) {
		//cout << it->first << "   ";  //debug
		tmp.push_back(it->first);
		it++;
	}
	return tmp;
}

/* return a vector of clusters as sets of members */
template<class T> void hatrees<T>::clusterArray(vector< set<T> > &vecset) {
	if (result.empty()) transform();
	vecset.clear();  // just in case it has something already in it

	typename std::multimap<T, T>::const_iterator mi = result.begin();
	while (mi != result.end()) {
		set<T> tmpset;
		tmpset.insert(mi->second);
		//os << mi->second;
		T cluster_id = mi->first;
		mi++;
		while (mi != result.end() && mi->first == cluster_id) {
			//os << '\t' << mi->second;
			tmpset.insert(mi->second);
			mi++;
		}
		vecset.push_back(tmpset);
		//os << endl;
	}
}

// you can provide a inital id to name the clusters
// This function will increment this id for each cluster
template<class T> void hatrees<T>::showClusterIntId(ostream &ous, int &id) {
	if (result.empty()) transform();
	typename std::multimap<T, T>::const_iterator mi = result.begin();
	while (mi != result.end()) {
		set<T> tmpset;
		ous << mi->second << '\t' << id << endl;
		//tmpset.insert(mi->second);
		T cluster_id = mi->first;
		mi++;
		while (mi != result.end() && mi->first == cluster_id) {
			//tmpset.insert(mi->second);
			ous << mi->second << '\t' << id << endl;
			mi++;
		}
		//vecset.push_back(tmpset);
		//os << endl;
		++id;
	}
}

/* output in cluster_id -->members 
 * Table format, good for relational database
 * */
template<class T> void hatrees<T>::showCluster(ostream &os, bool reverse) {
	if (result.empty()) transform();

	typename std::multimap<T, T>::const_iterator mi = result.begin();
	if (reverse) {
		while (mi != result.end()) {
			os << mi->second << '\t' << mi->first << endl;
			mi++;
		}
	}
	else {
		while (mi != result.end()) {
			os << mi->first << '\t' << mi->second << endl;
			mi++;
		}
	}
	//os << "Total members: " << result.size() << endl;
	//for db loading, this line causes trouble
	cerr << "Total members: " << result.size() << endl;
}

template<class T> void hatrees<T>::transform() {
	niterator it = nodes.begin();
	node<T>* root;
	while (it != nodes.end()) {
		root = getRoot(it->second->parent);
		result.insert(pair<T, T>(root->key, it->first));
		it++;
	}
}

/** Display one cluster per line, with all members on the same line
 * */
template<class T> int hatrees<T>::showClusterByLine(ostream &os) {
	if (result.empty()) transform();
	os << "Output format: one cluster per line\n\n";
	int clusterCnt=0;
	typename std::multimap<T, T>::const_iterator mi = result.begin();
	while (mi != result.end()) {
		clusterCnt++;
		os << mi->second;
		T cluster_id = mi->first;
		mi++;
		while (mi != result.end() && mi->first == cluster_id) {
			os << '\t' << mi->second;
			mi++;
		}
		os << "\n\n";
	}
	os << "Total " << clusterCnt << " clusters\n";
	return clusterCnt;
}

#endif
