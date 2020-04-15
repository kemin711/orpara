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

namespace orpara {
//using namespace pqxx;

/**
 * Node of the tree
 */
template<class T> struct HNode
{
	HNode() :  key(), parent(this), rank(0)  { }  // point to self
	HNode(const T &k) : key(k), parent(this), rank(0) { }
	HNode(T&& k) : key(std::move(k)), parent(this), rank(0) { }
	HNode(const T &k, HNode* p) : key(k), parent(p), rank(0) { }
	HNode(T&& k, HNode* p) : key(std::move(k)), parent(p), rank(0) { }
	/* set parent = this is essential */
	HNode(const HNode &n) : key(n.key), parent(this), rank(n.rank) { }
	HNode(HNode&& n) : key(std::move(n.key)), parent(this), rank(n.rank) { }
	// if use parent = n.parent then segmentation fault
	~HNode() { /*nothing needs to be done*/ }
	bool operator==(const HNode& n) const { return key==n.key; }
	bool operator<(const HNode &n) const { return key < n.key; }
	//HNode& operator=(const HNode &n) { if (&n != this) { key = n.key; parent = n.parent; rank = n.rank; } return *this; }
	HNode& operator=(const HNode &n) { if (&n != this) { key = n.key; parent = this; rank = n.rank; } return *this; }
	HNode& operator=(HNode&& n) { 
      if (&n != this) { 
         key = std::move(n.key); 
         parent = this; rank = n.rank; 
      } 
      return *this; 
   }

   /**
    * The look up key
    */
	T key;   // stored in the map<key, value>
   /**
    * Pointer to parent
    */
	HNode<T>* parent;
   /**
    * For algorithm operation
    */
	int rank;
};

/* tree operations */
///////////////////////////////////////////////////////
/* x is not altered */
template<class T> HNode<T>* getRoot(HNode<T> *x) {
	if (x != x->parent)  x->parent = getRoot(x->parent); 
	return x->parent;
}
template<class T> void mergeRoot(HNode<T> *&r1, HNode<T> *&r2) {
	if(r1 == r2) return;  // if both are the same set
	if (r1->rank > r2->rank) r2->parent = r1;
	else {
		r1->parent = r2;
		if (r1->rank == r2->rank) r2->rank++;
	}
}
template<class T> void join(HNode<T> &x, HNode<T> &y) {
	HNode<T> *r1 = getRoot(&x);  // this is necessary for the compiler
	HNode<T> *r2 = getRoot(&y);
	mergeRoot(r1, r2);
	//mergeRoot(getRoot(&x), getRoot(&y));  // compiler got confused 
}
template<class T> void join(HNode<T>* x, HNode<T>* y) {
	HNode<T> *r1 = getRoot(x);  // this is necessary for the compiler
	HNode<T> *r2 = getRoot(y);
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
//   Hatrees class definition /////////////////////////
/** 
 * this class provides simple file output methods.
 * No database table output is defined in this class
 * bedause it will require different database drivers.
 * The getCluster() method will produce a relational table.
 * The template type T must have the output operator<<
 * implemented to use the display functions.
 * For large object the T can ba pointer type.
 */
template<class T> class Hatrees
{
	public:
      /**
       * Default constructor
       */
		Hatrees() : nodes(), result() { }
      /** 
       * build a cluster object directly from multimap
       * Save the call to readFromMap()
       * @param mm is all the relationships among members.
       */
		Hatrees(const std::multimap<T,T> &mm)
         : nodes(), result()
      { readFromMap(mm); }
		~Hatrees();
		void clear() { nodes.clear(); result.clear(); }

		///////// input operations /////////////////
      /**
       * If you have an empty object then you can
       * take input from a multimap
       */
		void readFromMap(const std::multimap<T, T> &m);
      /**
       * Directly read from file for convenience.
       */
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

		/** 
       * changes the parent pointer of the nodes	
       */
		void showStore(ostream &os) const;  

		/** 
       * number of unique members 
       */
		int getNodeCount() const { return nodes.size(); }
		/** 
       * in a table format: (cluster_id(representative), members) 
       * */
		void showCluster(ostream &os, bool reverse= true);
		/** 
       * an initial id is given, the program will
       * auto_increment this id. This should work with
       * either string or int type of idtype.
       * Output member <TAB> id 
       * one row per member. 
       * */
		void showClusterIntId(ostream &ous, int &id);

		/** 
       * Display the result to output stream.
       * @return the number of clusters 
       * */
		int showClusterByLine(ostream &os);
      /** 
       * @return all the keys as a set<T> it is all the 
       *   members in the input. Useful for finding singleton
       *   clusters. 
       */
		set<T> keyset() const;
		vector<T> keyarray() const;

		/** 
       * produce an array of clusters, used call by reference
		 * to avoid creating tmp and copying.
       * The result could be used for further analysis. 
       * This is one of the output methods.
       * */
		void clusterArray(vector<set<T> > &vecset);

      /** 
       * @return a relational table
       *   with two columns: | representative | member |
       *   Use this method to load into database tables.
       *   Singletons are not included in the output.
       */
		//std::multimap<T,T>& getCluster() { 
		vector<pair<T,T>> getCluster() const { 
         if (result.empty()) transform();
         vector<T,T> res;
         for (auto& r : result) {
            for (auto& m : r.second) {
               res.push_back(make_pair(r.first->key, m));
            }
         }
         return res;
      }
      /**
       * @return a const reference to the result member.
       */
		const std::map<T*, vector<T>>& getResult() const { 
         if (result.empty()) transform(); 
         return result; 
      }
      int getNumberOfCluster() const {
         if (result.empty()) transform();
         return result.size();
      }
      vector<vector<T>> getClusterAsVector() const;
      vector<set<T>> getClusterAsSet() const;

#ifdef USE_HASH_MAP
      typedef typename hash_map<T, HNode<T>* >::iterator niterator;
      typedef typename hash_map<T, HNode<T>* >::const_iterator const_niterator;
#else
      typedef typename map<T, HNode<T>* >::iterator niterator;
      typedef typename map<T, HNode<T>* >::const_iterator const_niterator;
#endif

	private:
#ifdef USE_HASH_MAP 
		hash_map<T, HNode<T>* > nodes;
#else
		map<T, HNode<T>* > nodes;
#endif
		/** 
       * This is the result table:
       * representative->members 
       * Should have another method to return a reversed result.
       * Essentially the same as nodes, but more easy to work with.
       */
		//mutable map<T, vector<T>> result;
		mutable map<HNode<T>*, vector<T>> result;

		/** 
       * transform the nodStore cluster into mutimaped cluster
		 * in a table format cluster_id -> members 
       * Maybe more efficient to store member -> representative direction.
       * */
		void transform() const;
};

////////// Hatrees class template methods //////////////////////////

template<class T> Hatrees<T>::~Hatrees() { 
	niterator MI = nodes.begin();
	while (MI != nodes.end()) {
		delete MI->second;
		MI++;
	}
}

template<class T> void Hatrees<T>::readFromMap(const std::multimap<T, T> &m) {
	pair<niterator, bool> p1, p2;
   HNode<T>* tmpPtr=nullptr;
   auto imi = m.begin();
	while (imi != m.end()) {
      tmpPtr = new HNode<T>(imi->first);
		p1 = nodes.insert(make_pair(imi->first, tmpPtr));
      if (!p1.second) delete tmpPtr;
      tmpPtr = new HNode<T>(imi->second);
		p2 = nodes.insert(make_pair(imi->second, tmpPtr));
      if (!p2.second) delete tmpPtr;
		join(p1.first->second, p2.first->second);
		imi++;
	}
}

template<class T> void Hatrees<T>::readFromFile(const string &file) {
	ifstream IN(file.c_str());
	if (IN.fail()) {
		cerr << "reading from " << file << "failed\n";
		exit(1);
	}
	string ln;
	pair<niterator, bool> p1, p2;
	T f1, f2; 
   HNode<T>* tmpPtr;

	getline(IN, ln);
	while (!IN.eof()) {
		istringstream ist(ln);
		ist >> f1 >> f2;
      tmpPtr = new HNode<T>(f1);
		p1 = nodes.insert(make_pair(f1, tmpPtr));
      if (!p1.second) delete tmpPtr;
      tmpPtr = new HNode<T>(f2);
		p2 = nodes.insert(pair<T, HNode<T>* >(f2, tmpPtr));
      if (!p2.second) delete tmpPtr;
		join(p1.first->second, p2.first->second);
		getline(IN, ln);
	}
}

/* members->cluster_id sorted according to members if using <map>
 * if using hash_map then not in any order */
template<class T> void Hatrees<T>::showStore(ostream &os) const {
	const_niterator it = nodes.begin();
	HNode<T>* root;
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
template<class T> set<T> Hatrees<T>::keyset() const {
	set<T> tmpset;
	const_niterator it = nodes.begin();
	while (it != nodes.end()) {
		tmpset.insert(it->first);
		it++;
	}
	return tmpset;
}

template<class T> vector<T> Hatrees<T>::keyarray()  const {
	vector<T> tmp;
	const_niterator it = nodes.begin();
	while (it != nodes.end()) {
		tmp.push_back(it->first);
		it++;
	}
	return tmp;
}

template<class T> void Hatrees<T>::transform() const {
	const_niterator it = nodes.begin();
	HNode<T>* root;
	while (it != nodes.end()) {
		root = getRoot(it->second->parent);
		//result[root->key].push_back(it->first);
		result[root].push_back(it->first);
		it++;
	}
}

/* return a vector of clusters as sets of members */
template<class T> void Hatrees<T>::clusterArray(vector<set<T>> &vecset) {
	if (result.empty()) transform();
	vecset.clear();  // just in case it has something already in it
   for (auto& r : result) {
		set<T> tmpset(r.second.begin(), r.second.end());
		vecset.push_back(std::move(tmpset));
	}
}

template<class T> 
vector<vector<T>> Hatrees<T>::getClusterAsVector() const
{
   if (result.empty()) transform();
	vector<vector<T>> res;
   for (auto& r : result) {
      res.push_back(r.second);
	}
   return res;
}

template<class T> 
vector<set<T>> Hatrees<T>::getClusterAsSet() const {
	if (result.empty()) transform();
   vector<set<T>> res;
   for (auto& r : result) {
		res.push_back(set<T>(r.second.begin(), r.second.end()));
	}
   return res;
}


// you can provide a inital id to name the clusters
// This function will increment this id for each cluster
template<class T> void Hatrees<T>::showClusterIntId(ostream &ous, int &id) {
	if (result.empty()) transform();
   for (auto& r : result) {
      for (auto& x : r.second) {
         ous << x << '\t' << id << endl;
      }
		++id;
	}
}

/* output in cluster_id -->members 
 * Table format, good for relational database
 * */
template<class T> void Hatrees<T>::showCluster(ostream &os, bool reverse) {
	if (result.empty()) transform();
   int memberCnt=0;
	if (reverse) {
      for (auto& r : result) {
         for (auto& m : r.second) {
            os << m << '\t' << r.first->key << endl;
            ++memberCnt;
         }
		}
	}
	else {
      for (auto& r : result) {
         for (auto& m : r.second) {
            os << r.first->key << '\t' << m << endl;
            ++memberCnt;
         }
		}
	}
	cerr << "Total members: " << memberCnt << endl;
}

/* Display one cluster per line, with all members on the same line
 * */
template<class T> int Hatrees<T>::showClusterByLine(ostream &os) {
	if (result.empty()) transform();
	os << "Output format: one cluster per line\n\n";
   for (auto& r : result) {
      os << r.second[0];
      for (int i=1; i<r.second.size(); ++i) {
         os << '\t' << r.second[i];
      }
      os << endl;
   }
	os << "Total " << getNumberOfCluster() << " clusters\n";
	return getNumberOfCluster();
}

} // end of orpara namespace
#endif
