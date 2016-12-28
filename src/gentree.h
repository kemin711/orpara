#ifndef GENTREE_H
#define GENTREE_H

//#include <iostream>
//#include <fstream>
//#include <string>
#include <queue>

using namespace std;

namespace orpara {
// count the leading space
// int leadingspace(const string &s);
// this is a header only template header there is no need for linking
/** file: gentree.h
 * class node
 * */
template<class T> struct node {
	/**  
    * Default constructor.
    * Parent is self, child and sibling are 0.
    */
	node() : data(), child(0), sibling(0), parent(this) {}

	node(const T& val) : child(0), sibling(0), parent(this), data(val) { }
	/** 
    * Constructor from data value and parent
    * @param val the data value of the node
    * @param p the parent pointer
    */
	node(const T& val, node<T> *p) : child(0), sibling(0), data(val), parent(p) { }
   // not sure you need to have a destructor
   //~node();

	// copy constructor does not make sense for tree
	//node(const node<T>& n) { data=n.data; child=n.child; sibling=n.sibling; parent=n.parent; } 

	/** 
    * only the value of the node is needed
    * Use (value, parent) constructor
    */
	void add_child(const T& chi) { child = new node<T>(chi, this); }
	void addChild(const T& chi) { child = new node<T>(chi, this); }
   /**
    * apply this function to this node 
    * and the last sibling of this node to build a chain of siblings.
    */
	void add_sibling(const T& sib) { sibling = new node<T>(sib, this->parent); }
   /**
    * If no sibling yet add sibling.
    * If there are siblings, then append to the tail of
    * the sibling chain.
    */
   node<T>* addSibling(const T& sib);
   /**
    * Add to the end of the chain of siblings
    */
   void append_sibling(const T& sib);
	void setData(const T& d) { data = d; }
	node<T>* right() { return sibling; }
	node<T>* down() { return child; }
   /**
    * Get all the siblings
    */
	vector<T> row(); 
   /**
    * Get all the first child
    */
	vector<T> col();
   bool isRoot() const { return parent == this; }

	T data;   // data stored in this node
	node<T> *child;   // left or down (in ACE) child
	node<T> *sibling; // righ sibling
	node<T> *parent;
};

//typedef template<class T> void (*visitFunc)(node<T> *n);

template<class T> class GenericTree {
   public:
      /**
       * Default constructor build an empty tree with root node only.
       */
      GenericTree() : root(new node<T>()) {}
      void deallocate();
      ~GenericTree() { deallocate(); }
      node<T>* root;
};

////////////// not class implementation ///////////////////////

template<class T> node<T>* node<T>::addSibling(const T& sib) { 
   if (sibling == 0) {
      sibling = new node<T>(sib, this->parent);
      return sibling;
   }
   else {
      node<T>* p = sibling;
      while (p->sibling != 0) p=p->sibling;
      p->sibling = new node<T>(sib, this->parent);
      return p->sibling;
   }
}

template<class T> void node<T>::append_sibling(const T& sib) {
   node<T> p=this;
   while (p->sibling != 0) p = p->sibling;
   p->sibling = new node<T>(sib, this->parent);
}

template<class T> vector<T> node<T>::row() { 
	node<T>* ptr = sibling;
	vector<T> tmp;
	while (ptr) { 
		tmp.push_back(ptr->data);
		ptr = ptr->sibling;
	}
	return tmp;
}

template <class T> vector<T> node<T>::col() {
	node<T>* ptr = sibling;
	vector<T> tmp;
	while (ptr) {
		tmp.push_back(ptr->data);
		ptr = ptr->child;
	}
	return tmp;
}

//template<class T> void preorder(node<T> *n, visitFunc vis) {
template<class T> void preorder(node<T> *n, void (*vis)(node<T>*)) {
	if (n) {
		vis(n);
		//cout << n->data << " || ";
		preorder(n->child, vis);
		preorder(n->sibling, vis);
	}
	else {
		vis(0); //cout << endl;
	}
}

//template<class T> void levelorder(node<T> *n, visitFunc vis) {
/* visit the horizontal siblings first
 */
template<class T> void levelorder(node<T> *n, void (*vis)(node<T>*)) {
	queue<node<T>* > nodes_left;
	nodes_left.push(n);
	node<T>* p;

	while (!nodes_left.empty()) {
		p = nodes_left.front();
		//vis(p->parent);   //cout << "parent:" << p->parent->data << endl;
		while (p) {
			vis(p);  //cout << p->data << " || ";
			if (p->child) nodes_left.push(p->child);
			p = p->sibling;
		}
		vis(0);  //cout << endl << endl;
		nodes_left.pop();
	}
}

//template<class T> void postorder(node<T> *n, visitFunc vis) {
template<class T> void postorder(node<T> *n, void (*vis)(node<T>*)) {
	if (n) {
		postorder(n->child, vis);
		postorder(n->sibling, vis);
		vis(n);  //cout << n->data << " || ";
	}
	else {
		vis(0);  //cout << endl;
	}
}

// postorder traversal to delete, the only way of 
// deleting all nodes
template<class T> void deltree(node<T> *n) { 
   if (n) {
       deltree(n->child);
       deltree(n->sibling);
       //cout << "deleting " << n->data << endl;  //debug
       delete n;
   }
}

///////////////// GenericTree class /////////////
template<class T> void GenericTree<T>::deallocate() { 
   deltree(root);
}
}

#endif
