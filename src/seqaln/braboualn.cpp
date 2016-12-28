#include "braboualn.h"
#include <iostream>
//#include <array>
#include <stack>

namespace orpara {
Node::Node(Node* pn, AlignType at) : p1(pn->p1), p2(pn->p2),
      align(at), PrS(pn->PrS), Tmin(pn->Tmin), Tmax(pn->Tmax),
      parent(pn), child{0,0}, depth(pn->depth+1),
      identical(pn->identical)
{ 
   if (align == MATCH) {
      ++identical;
      ++p1; ++p2;
   }
   else if (align == MISMATCH) {
      ++p1; ++p2;
   }
   else if (align == S1GAP) ++p2;
   else if (align == S2GAP) ++p1;
   else { // for debug
      cerr << "wrong align type " << at << endl;
      exit(1);
   }
}

void Node::setChildren(Node* n1, Node* n2) {
   if (child[0] != 0 || child[1] != 0) {
      cerr << "parent node already has children, review algorithm!\n";
      exit(1);
   }
   child[0]=n1; child[1]=n2; 
}

void Node::setChildren(Node* n1, Node* n2, Node* n3) {
   if (child[0] != 0 || child[1] != 0 || child[2] != 0) {
      cerr << "parent node already has children, review algorithm!\n";
      exit(1);
   }
   child[0]=n1; child[1]=n2; child[2]=n3; 
}

void Node::setChildren(const array<Node*, 3> &nodes) {
   if (child[0] != 0 || child[1] != 0 || child[2] != 0) {
      cerr << "parent node already has children, review algorithm!\n";
      exit(1);
   }
   child[0]=nodes[2];
   child[1]=nodes[1];
   child[2]=nodes[0];
}

ostream& operator<<(ostream &ous, const Node &n) {
   ous << "p1: " << n.p1 << " p2 " << n.p2 << endl
      << "PrS: " << n.PrS << " Tmin: " << n.Tmin << " Tmax: " << n.Tmax
      << endl;
   switch(n.align) {
      case MATCH:
         ous << "match";
         break;
      case MISMATCH:
         ous << "mismatch";
         break;
      case S1GAP:
         ous << "seq1 gap";
         break;
      case S2GAP:
         ous << "seq2 gap";
         break;
      default:
         ous << "default match type";
   }
   ous << endl;
   ous << "identity: " << n.getIdentity() << endl;
   if (n.parent == 0) ous << "Root node\n";
   return ous;
}

////// fitness queue ////////
FitnessQueue::FitnessQueue(int seq1len, int seq2len, int match, int mism, int gap)
   : basei(min(seq1len,seq2len)*mism + abs(seq1len-seq2len)*gap),
      items(min(seq1len,seq2len)*match+abs(seq1len-seq2len)*gap-basei+1),
      maxi(0)
{
   cout << "number of items: " << items.size() << endl;
}

FitnessQueue::FitnessQueue(int seq1len, int seq2len, int match, int mism, int gapo, int gape)
   : basei(min(seq1len,seq2len)*mism + gapo*abs(seq1len-seq2len)),
      items(min(seq1len,seq2len)*match+gapo+gape*(abs(seq1len-seq2len)-1)-basei+1),
      maxi(0)
{
   //cout << "Priority queue capacity: " << items.size() << endl;
}

void FitnessQueue::enqueue(Node* node) {
   int i = node->Tmax - basei;
   //cout << "queue location " << i << endl;
   items[i].insert(node);
   if (i>maxi) maxi=i;
}

Node* FitnessQueue::dequeue() {
   Node * topNode = *(items[maxi].begin());
   items[maxi].erase(items[maxi].begin());
   while (items[maxi].empty()) {
      --maxi;
   }
   return topNode;
}

void FitnessQueue::show() const {
   cout << "basei: " << basei
      << " maxi: " << maxi
      << " table size: " << items.size() << endl;
}

////// BraBoualn //////////////
pair<int, int> Braboualn::futureScore(int x1, int x2) {
   int shorter=min(x1,x2);
   int diff=abs(x1-x2);
   return make_pair(shorter*Ss + Go*diff, shorter*Sm + affine(diff));
}

//      int matchS=5, int mismatchS=-4, int gapO=-9, int gapE=-1) 
Braboualn::Braboualn(const string &s1, const string &s2, int matchS, int mismatchS, int gapO, int gapE) 
      : seq1(s1), seq2(s2), cutoff(0.3), Sm(matchS), Ss(mismatchS), 
      Go(gapO), Ge(gapE),
      candidate(s1.length(), s2.length(), Sm, Ss, Go, Ge), 
      root(new Node()),
      best(0), current(root), explored(s1.length()+1, vector<Node*>(s2.length()+1,0)),
      alignedseq()
{ 
   root->setFitness(futureScore(s1.length(), s2.length()));
   //cout << "after constructor:\n";
   //show();
}

void deallocateNode(Node *n) {
   if (n==0) return;
   deallocateNode(n->child[0]);
   deallocateNode(n->child[1]);
   deallocateNode(n->child[2]);
   delete n;
   n=0;
}

Braboualn::~Braboualn() {
   deallocateNode(root);
}

// produce the three children in order of Tmax from large to small
// expand the current node
// The boundary condition needs to be considered
// dealing with the last (or leaf) residues.

int Braboualn::affine(int g) const {
   if (g<0) {
      cerr << "warnning gap value less than zero: " << g << endl;
   }
   if (g <= 0) return 0;
   if (g == 1) return Go;
   return Go + (g-1)*Ge;
}

Node* Braboualn::produceGap1Child() {
   // these check should be moved outside of a main loop
   // for performance. For the first version we are fine.
   if (current->align == S2GAP) {
      //cerr << "cannot switch from GAP in seq2 to gap in seq1\n";
      return 0;
   }
   //make sure seq2 is not at last residue
   if (current->p2 == seq2.length()-1) { // cannot extend gap1
      if (current->align == MATCH || current->align == MISMATCH) {
         //cerr << "seq2 at last char, cannot produce gap1 child\n";
         return 0;
      }
   }
   Node* gap1Child = new Node(current, S1GAP);
   if (current->align == MATCH || current->align == MISMATCH || current->parent == 0) {
      gap1Child->PrS += Go;
   }
   else if (current->align == S1GAP) {
      gap1Child->PrS += Ge;
   }
   else { // exception, should not go here, remove in production
      cerr << "align type of current unknown: " << current->align << endl;
      exit(1);
   }
   int R1=seq1.length()-current->p1;
   int R2=seq2.length()-current->p2;
   int Fmin, Fmax;
   Fmin=min(R1, R2-1)*Ss + Go*abs(R1-R2+1);
   Fmax=min(R1, R2-1)*Sm + affine(abs(R1-R2+1));
   gap1Child->Tmin=gap1Child->PrS + Fmin;
   gap1Child->Tmax=gap1Child->PrS + Fmax;
   return gap1Child;
}

Node* Braboualn::produceGap2Child() {
   if (current->align == S1GAP) {
      //cerr << "cannot switch from GAP in seq1 to gap in seq2\n";
      return 0;
   }
   // seq1 is not at last residue
   if (current->p1 >= seq1.length()) {
      //cerr << "seq at last char, cannot produce gap2 child\n";
      return 0;
   }
   Node* gap2Child = new Node(current, S2GAP);
   if (current->align == MATCH || current->align == MISMATCH || current->parent == 0) {
      gap2Child->PrS += Go;
   }
   else if (current->align == S2GAP) {
      gap2Child->PrS += Ge;
   }
   else {
      cerr << "current bad align type\n";
      exit(1);
   }
   int R1=seq1.length()-current->p1;
   int R2=seq2.length()-current->p2;
   int Fmin, Fmax;
   Fmin=min(R1-1, R2)*Ss + Go*abs(R1-1-R2);
   Fmax=min(R1-1, R2)*Sm + affine(abs(R1-1-R2));
   gap2Child->Tmin=gap2Child->PrS + Fmin;
   gap2Child->Tmax=gap2Child->PrS + Fmax;
   return gap2Child;
}

Node* Braboualn::produceAlignChild() {
   /* checked by the two other methods.
   if (current->p1 == seq1.length()-1 || current->p2 == seq2.length()-1) {
      cerr << "there is no child for terminal current node\n";
      return 0;
   }
   */
   Node *alignChild;
   //cout << "seq1,2 char: " << P1 << ":" << seq1[P1] << "x" 
   //   << P2 << ":" << seq2[P2] << endl;
   if (seq1[current->p1] == seq2[current->p2]) { // match case
      alignChild = new Node(current, MATCH);
      alignChild->PrS += Sm;
      alignChild->Tmin += (Sm-Ss);
   }
   else {
      alignChild = new Node(current, MISMATCH);
      alignChild->PrS += Ss;
      alignChild->Tmax += Ss-Sm;
   }
   return alignChild;
}

void Braboualn::giveBirth() {
   //cerr << "giveBirth()\n";
   // P1 and P2 should not exceed the last index
   int P1=current->p1;
   int P2=current->p2;
   /*cout << "P1, P2: " << P1 << ", " << P2 << endl;
   if (P1 < seq1.length() && P2 < seq2.length()) {
      cout << "seq 1, 2 char: " << seq1[P1] << ", " << seq2[P2] << endl;
   }
   else {
      cout << "at end of one of the two sequences\n";
      if (P1 < seq1.length()) {
         cout << "seq1 char: " << seq1[P1] << endl;
      }
      if (P2 < seq2.length()) {
         cout << "seq2 char: " << seq2[P2] << endl;
      }
   }
   cout << "last index of seq: " << seq1.length()-1
      << ", " << seq2.length()-1 << endl;
      */
   // this is the leaf node special case
   // Boundary condition
   // seq2 reach the end, seq1 still has residues
   if (P2 >= seq2.length()) {
      if (P1 < seq1.length()) {
         current->setFirstChild(produceGap2Child());
         return;
      }
      else {
         cerr << "both seq1 and seq2 reach the ends\n";
         exit(1);
      }
   }
   if (P1 >= seq1.length()) {
      if (P2 < seq2.length()) {
         current->setFirstChild(produceGap1Child());
         return;
      }
      else {
         cerr << "both seq1 and seq2 reach the ends\n";
         exit(1);
      }
   }

   /*
   if (P1 >= seq1.length()-2 && P2 >= seq2.length()-2) {
      current->setFirstChild(produceAlignChild());
      return;
   }
   if (P1 >= seq1.length()-1) { // seq1 reached the last residue
      // you can only produced gapped alignment after this
      current->setFirstChild(produceGap1Child());
      return;
   }
   if (P2 >= seq1.length()-1) { // seq2 last residue
      current->setFirstChild(produceGap2Child());
      return;
   }
   */
   // child with match, or gap in the short sequence
   // always have the largest Tmax
   Node *alignChild=produceAlignChild();
   Node *gap1Child=produceGap1Child();
   Node *gap2Child=produceGap2Child();
   if (gap1Child == 0 && gap2Child == 0) {
      current->setFirstChild(alignChild);
      return;
   }
   if (gap1Child == 0) {
      if (alignChild->betterThan(gap2Child)) {
         current->setChildren(alignChild, gap2Child);
      }
      else {
         current->setChildren(gap2Child, alignChild);
      }
      return;
   }
   if (gap2Child == 0) {
      if (alignChild->betterThan(gap1Child)) {
         current->setChildren(alignChild, gap1Child);
      }
      else {
         current->setChildren(gap1Child, alignChild);
      }
      return;
   }
   //cout << "gap2Child: " << *gap2Child << endl;
   array<Node*, 3> tmp = {alignChild, gap1Child, gap2Child};
   sort(tmp.begin(), tmp.end(), sortbyTmax());
   /*
   cout << "after sorting of three children:\n"
      << "-----------------\n"
      << *tmp[0] << "-----------------\n"
      << *tmp[1] << "-----------------\n"
      << *tmp[2] << "-----------------\n" << endl;
      */
   current->setChildren(tmp[2], tmp[1], tmp[0]);
}

/*
double Braboualn::getImax() const {
   double bestIdentity = (double)(current->identical + getShorterTailLength())/(current->depth + getLongerTailLength()); 
   //cout << "best identity achievable: " << bestIdentity << endl;
   return bestIdentity;
}
*/

bool Braboualn::promising(const Node* n) const {
   if (isLeafNode(n)) {
      if (best == 0) return true;
      if (n->Tmax > best->Tmax) return true;
      return false;
   }
   if (explored[n->p1][n->p2] == 0) { // no same node to compare
      if (best == 0) return true;
      if (n->Tmax > best->Tmax) return true;
      return false;
   }
   //cout << explored[n->p1][n->p2]->p2 << endl;
   //cout << "explored PrS " << explored[n->p1][n->p2]->PrS << endl;
   if (n->PrS <= explored[n->p1][n->p2]->PrS) { return false; }
   if (best == 0) return true;
   if (n->Tmax > best->Tmax) return true;
   return false;
}

void Braboualn::enqueueSibling() {
#ifdef IGNORE_THIRD_CHILD
   for (int i=1; i<2; ++i) {
#else
   for (int i=1; i<3; ++i) {
#endif
      if (current->getChild(i) == 0) continue;
      if (promising(current->getChild(i))) {
         candidate.enqueue(current->getChild(i));
      }
      else { current->pruneChild(i); }
   }
}

void Braboualn::showExplored() const {
   //cout << string(30, '=') << endl;
   //cout << "Explored nodes so far:\n";
   for (size_t i=0; i<explored.size(); ++i) {
      int numitem=0;
      for (size_t j=0; j<explored[i].size(); ++j) {
         if (explored[i][j] != 0) {
            cout << " " << *explored[i][j];
            ++numitem;
         }
      }
      if (numitem > 0) cout << string(40, '-') << endl;
   }
}

int Braboualn::run() {
   int P1, P2;
   int count=0; // for algorithm investigation
   // this will be removed after the release stage.
   do {
      P1=current->p1; P2=current->p2;
      // not leaf node yet
      while (P1 < seq1.length() || P2 < seq2.length()) {
         ++count;
         //cout << "\ncurrent node: " << *current << endl;
         // expand until node at end of at least one of the sequences
         giveBirth();
         enqueueSibling();
         Node* C1 = current->getFirstChild();
         //cout << C1 << endl;
         //cout << "child 1\n" << *C1 << endl;
         explored[P1][P2]=current;
         if (!promising(C1)) {
            current->pruneFirstChild();
            break;
         }
         current=C1;
         P1=C1->p1; P2=C1->p2;
      }
      if (P1 == seq1.length() && P2 == seq2.length()) { // child is leaf now
         if (best == 0) best=current;
         else if (current->Tmax > best->Tmax) {
            best = current;
         }
         //cout << "best so far:\n" << *best;
      }
      current = candidate.dequeue();
      if (getImax() < cutoff) { // cutoff 30% by default
         //cout << "Maximum identity lower than " << cutoff 
         //   << ", giving up current:\n"
         //   << *current << endl;
         break;
      }
   } while (current->Tmax >= best->Tmax);
   //cout << "best score: " << *best << endl;
   traceBack(best);
   //cout << "total nodes examined: " << count << endl;
   return count;
}

void Braboualn::show() const {
   cout << "seq1 and 2:\n" 
      << seq1.length() << " " << seq2.length() << endl
      << seq1 << endl << seq2 <<endl
      << "Parameters Sm " << Sm << " Ss " << Ss << " Gap " << Go << ":" << Ge << endl
      << "root: ";
   if (root == 0) { cout << 0 << endl; }
   else { cout << *root; }
   cout << "best: ";
   if (best == 0) cout << 0 << endl;
   else cout << *best;
   cout << "current: ";
   if (current == 0) cout << 0 << endl;
   else cout << *current;
   cout << "candidate: ";
   candidate.show();
   cout << "========" << endl;
}

void Braboualn::traceBack(Node* n) {
   string top(n->depth, '-'), bottom(n->depth, '-');
   size_t i=n->depth-1;
   while (n->parent != 0) {
      //cout << string(30, '-') << endl;
      //cout << *n << endl << string(30, '-') << endl;
      if (n->align == MATCH || n->align == MISMATCH) {
         top[i] = seq1[n->p1-1];
         bottom[i] = seq2[n->p2-1];
      }
      else if (n->align == S1GAP) {
         bottom[i] = seq2[n->p2-1];
      }
      else if (n->align == S2GAP) {
         top[i] = seq1[n->p1-1];
      }
      else {
         cerr << "match type: " << n->align << " match type unknown\n";
         cerr << "must be a legal match type\n";
         exit(1);
      }
      n = n->parent;
      --i;
   }
   //cout << "alignment:" << top.size() << " alnlen\n";
   alignedseq=make_pair(std::move(top), std::move(bottom));
}

void Braboualn::displayAlignment(ostream& ous) const {
   size_t width=70;
   size_t i=0;
   ous << "identity=" << getIdentity() 
      << " alignlen=" << getAlnlen() << endl << endl;
   while (i < alignedseq.first.length()) {
      string str1(alignedseq.first.substr(i,width));
      string str2(alignedseq.second.substr(i,width));
      ous << str1 << endl;
      for (size_t j=0; j<str1.length(); ++j) {
         if (str1[j] == str2[j]) ous << '|';
         else ous << ' ';
      }
      ous << endl << str2 << endl << endl;
      i += width;
   }
}
}
