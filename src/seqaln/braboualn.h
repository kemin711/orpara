#ifndef BRANCHBOUND_H
#define BRANCHBOUND_H
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <iostream>
#include <array>
//if explore the third child, it produce slightly better
//alignments, but at the cost of 50% drop in performance.
//#define IGNORE_THIRD_CHILD

using namespace std;

namespace orpara {
// UNDEF is for the ROOT Node because there is no
// input yet.
enum AlignType { MATCH, MISMATCH, S1GAP, S2GAP, UNDEF };

class Node {
   public:
      Node() : p1(0), p2(0), align(UNDEF), PrS(0), Tmin(0), Tmax(0), 
            parent(0), child{0,0,0}, depth(0), identical(0) { }
      void setFitness(pair<int,int> range) { Tmin=range.first; Tmax=range.second; }
      /**
       * @apram pn is the parent node.
       * Construct a node from another node pointer
       * Only the depth got incremented. All other 
       * parameters are copied.
       */
      Node(Node* pn, AlignType at);
      //: p1(pn->p1), p2(pn->p2),
      //   align(at), PrS(pn->PrS), Tmin(pn->Tmin), Tmax(pn->Tmax),
      //      parent(pn), child({0,0}), depth(pn->depth+1),
      //      identical(pn->identical)
      //      { }
      bool isAlign() const { return align == MATCH || align == MISMATCH; }
      bool isIndel1() const { return align == S1GAP; }
      bool isIndel2() const { return align == S2GAP; }
      void setFirstChild(Node* node) { child[0] = node; }
      void setSecondChild(Node* node) { child[1] = node; }
      void setThirdChild(Node* node) { child[2] = node; }
      /**
       * Make sure there is no existing children!
       * Algorithm needs to be refined. Make sure no 
       * memory leaks.
       */
      void setChildren(Node* n1, Node* n2, Node* n3);
      void setChildren(Node* n1, Node* n2);
      void setChildren(const array<Node*,3> &nodes);
      //{ child[0]=n1; child[1]=n2; }
      Node* getSecondChild() const { return child[1]; }
      Node* getFirstChild() const { return child[0]; }
      Node* getChild(int i) const { return child[i]; }
      void pruneFirstChild() { delete child[0]; child[0]=0; }
      void pruneSecondChild() { delete child[1]; child[1]=0; }
      void pruneChild(int i) { delete child[i]; child[i]=0; }
      double getIdentity() const { if (depth==0) return 0; return (double)identical/depth; }
      friend ostream& operator<<(ostream &ous, const Node &n);
      bool betterThan(const Node* n) const {
         if (Tmax > n->Tmax) return true;
         if (Tmax < n->Tmax) return false;
         return Tmin > n->Tmin;
      }

   public:
      /** the index in the explored node containner.
       * The index in the input sequence is p1-1, or p2-1
       * p1 for the first sequence and p2 for the second sequence
       */
      int p1, p2;
      /** the alignment type, this may not be needed
       * p1 or p2 set to -1 may be used for a gap. */
      AlignType align; // type aligned or indel
      /** PrS, present score */
      int PrS;
      /** store the Tmin and Tmax */
      //pair<int, int> fitness;
      int Tmin, Tmax;
      /** for trace back */
      Node* parent;
      /** Child 1 for match, mismatch
       * child 2 for gap in S2
       * child 3 for gap in S1
       * Only expand maximum of two child
       */
      Node* child[3];
      /** distance to root, also the alignment length */
      int depth;
      /** number of identical residues on the path
       * for reading the identity of the alignment
       */
      int identical;
};

class sortbyTmin {
   public:
      bool operator()(const Node* node1, const Node* node2) {
         return node1->Tmin < node2->Tmin;
      }
};
class sortbyTminDesc {
   public:
      bool operator()(const Node* node1, const Node* node2) {
         return node1->Tmin > node2->Tmin;
      }
};

class sortbyTmax {
   public:
      bool operator()(const Node* node1, const Node* node2) {
         if (node1->Tmax < node2->Tmax) return true;
         if (node1->Tmax > node2->Tmax) return false;
         return (node1->Tmin < node2->Tmin);
      }
};

class FitnessQueue {
   public:
      /**
       * Constructor for the simple Gap score.
       */
      FitnessQueue(int seq1len, int seq2len, int match, int mism, int gap);
      /**
       * for affine gap score
       */
      FitnessQueue(int seq1len, int seq2len, int match, int mism, int gapo, int gape);
      //FitnessQueue(int bottom, int top) 
      //   : basei(bottom), items(top-bottom), maxi(0) {
      //   cout << "bottom: " << bottom << endl 
      //   << "top: " << top << endl; 
      //}
      void enqueue(Node* node);
      Node* dequeue();
      /** for debug */
      void show() const;

   private:
      int basei;
      /** remember the largest index for Tmax */
      vector< set<Node*, sortbyTminDesc> > items;
      int maxi;
};


/**
 * This class will be used an algorithm object to 
 * process multiple inputs for efficiency.
 */
class Braboualn {
   public:
      /**
       * Start from two sequences
       */
      //Braboualn(const string &s1, const string &s2);
      Braboualn(const string &s1, const string &s2,
                  int matchS=5, int mismatchS=-4, int gapO=-9, int gapE=-1);

      /** 
       * set the best pointer.
       * No need for return.
       * @return the number of nodes explored.
       * This is for algorithm investigation. Will
       * change the return type in beta version.
       */
      int run();
      /** produce three children from the current node
       * and order them in the array according to their 
       * Tmax score from small to large
       * Applicable only if current node
       * for both sequences are not at the last residue:
       * P1<seq1.length()-1 and P2<seq1.length()-1
       * Set the first two children of the current node
       * to the first two children according the Tmax.
       *
       * @return 3 child nodes from small to large
       */
      //array<Node*, 3> giveBirth();
      void giveBirth();
      /** the best identity the current node can achieve
       * It is the implied Tmax.
       */
      int getShorterLength() const { return min(seq1.length(), seq2.length()); }
      int getLongerLength() const { return max(seq1.length(), seq2.length()); }
      int getShorterTailLength() const { return min(seq1.length()-current->p1, seq2.length()-current->p2); }
      int getLongerTailLength() const { return max(seq1.length()-current->p1, seq2.length()-current->p2); }
      double getImax() const {
          return (double)(current->identical + getShorterTailLength())/(current->depth + getLongerTailLength()); }
      ~Braboualn();
      /** debug function */
      void show() const;
      /**
       * given a node, it produce the sequence alignment
       */
      void traceBack(Node* n);
      bool promising(const Node* n) const;
      void enqueueSibling();
      void showExplored() const;
      int affine(int g) const;
      void setGap(int gopen, int gext) { Go=gopen; Ge=gext; }
      void setIdentityCutoff(double idcut) { cutoff=idcut; }
      double getIdentity() const { if (best==0) return 0; return best->getIdentity(); }
      int getGapOpen() const { return Go; }
      int getGapExtend() const { return Ge; }
      void displayAlignment(ostream& ous) const;
      int getMatch() const { return Sm; }
      int getMismatch() const { return Ss; }
      int getAlnlen() const { return best->depth; }

   private:
      Node* produceAlignChild();
      Node* produceGap1Child();
      Node* produceGap2Child();
      bool isLeafNode(const Node* n) const { return n->p1 == seq1.length() && n->p2 == seq2.length(); }
      /** compute the future score given the 
       * length x1 and x2 for the remaining length of the 
       * not aligned sequences.
       * @return (Fmin, Fmax)
       */
      pair<int, int> futureScore(int x1, int x2);
      pair<int, int> futureScoreAffine(int x1, int x2);
      /** 
       * Future version should use reference to objects 
       * to reduced the copying cost
       * Use string to save the trouble of managing memory.
       */
      string seq1, seq2;
      /**
       * Identity cutoff, default 0.3
       */
      float cutoff;
      /** use simple number for the first version
       * fugure versions will use templated scoreing matrix
       */
      int Sm;
      /** mismatch score for one paire of residues */
      int Ss;
      /**
       * For affine gap score. This one gets better results
       * than the simple one.
       */
      int Go, Ge;
      /** simulate priority queue behavior */
      FitnessQueue candidate;
      /** the root node, this will store the 
       * hashing base index value too
       * Tmin is this value. */
      Node *root;
      /** pointer to the leaf of the best branch
       */
      Node *best;
      /** for algoritm execution */
      Node *current;
      /** all the nodes expanded for the whole algorithm
       * This server as a look up by the index of the two 
       * sequences.
       * for fast indexing of the nodes in the tree.
       * For performance, this can be replaced with 
       * two number look up hash
       * map<int, map<int>> explored;
       */
      vector<vector<Node*> > explored;
      pair<string,string> alignedseq;
};
}
#endif
