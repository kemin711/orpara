#include "dynaln.h"
#include <iomanip>
#include <sstream>
#include "strformat.h"

//#define DEBUG
// replaced by the template version all apps using this 
// should be updated

namespace orapra {
// helper function
// this is not needed
/*
int endgaps(const string &gs) {
   int i=gs.length()-1;
   int L=0;
   while (i>=0 && gs[i--] == Dynaln::gapchar) {
      ++L;
   }
   if (i==0) {
      cerr << "problem, the whole string is gap char!\n";
      exit(1);
   }
   return L;
}
*/

Dynaln::~Dynaln() {
   if (S != 0) delete[] S; 
   delete[] M;
   delete[] IX;
   delete[] IY;
}

void Dynaln::setGapParameterFromMatrix() {
   setGapInsert(ST.getGapInsert());
   setGapExtend(ST.getGapExtend());
}

//void Dynaln::allocmem(const int Nr, const int Nc) {
void Dynaln::allocmem() {
   int Nr=seq1->length();
   int Nc=seq2->length();
   if (S == 0) { // allocate from new memory
      Ssize=Nr*Nc;
      numcol=Nc;
      S = new int[Ssize];
      M = new int[Nc+1]; // for match state
      IX = new int[Nc+1]; // for Insertion of X
      IY = new int[Nc+1]; // for insertion of Y
   }
   else {
      if (Ssize < Nr*Nc) { // old memory not enough
         delete[] S;
         Ssize=Nr*Nc;
         S=new int[Ssize];
      }
      if (numcol < Nc) { // column memory not enough
         delete[] M; 
         delete[] IX; 
         delete[] IY;
         M = new int[Nc+1];
         IX = new int[Nc+1];
         IY = new int[Nc+1];
         numcol=Nc;
      }
   }
   // asigne the code of each sequences
   C1=seq1->getcode();
   C2=seq2->getcode();
}

/*
 * use last direction as a guide to break a tie
 * 0---j--->
 * |
 * i
 * |
 *\|/
 */
void Dynaln::tracepointer(int &i, int &j) {
   if (!alnidx.empty()) {
      alnidx.clear();
   }
   int Nc=seq2->length();
   int pi, pj;
   int *ptr = S+(i*seq2->length() + j);
   int lastTrace = PTRNULL; // should have only one direction
   while (i>-1 && j>-1) {
      if ( *ptr == PTRDIAG) {
         alnidx.push_front(make_pair(i,j));
         lastTrace = PTRDIAG;
         ptr -= (Nc+1);
         --i; --j;
      }
      else if (*ptr == PTRLEFT) {
         alnidx.push_front(make_pair(-1, j));
         lastTrace = PTRLEFT;
         --ptr;
         --j;
      }
      else if (*ptr == PTRTOP) {
         alnidx.push_front(make_pair(i, -1));
         lastTrace = PTRTOP;
         ptr -= Nc;
         --i;
      }
      else if ( ((*ptr)&PTRDIAG) == PTRDIAG && ((*ptr)&PTRLEFT) == PTRLEFT
            && ((*ptr)&PTRTOP) == PTRTOP) { // all three directions
         if (lastTrace == PTRDIAG) {
            alnidx.push_front(make_pair(i,j));
            ptr -= (Nc+1);
            --i; --j;
         }
         else if (lastTrace == PTRLEFT) {
            alnidx.push_front(make_pair(-1, j));
            --ptr;
            --j;
         }
         else { // PTRTOP
            alnidx.push_front(make_pair(i, -1));
            ptr -= Nc;
            --i;
         }
      }
      else if ( ((*ptr)&PTRDIAG) == PTRDIAG && ((*ptr)&PTRTOP) == PTRTOP) { // both top and diag
         if (lastTrace == PTRDIAG) {
            alnidx.push_front(make_pair(i,j));
            ptr -= (Nc+1);
            --i; --j;
         }
         else if (lastTrace == PTRTOP) {
            alnidx.push_front(make_pair(i, -1));
            ptr -= Nc;
            --i;
         }
         else { // just pick diagnal
            alnidx.push_front(make_pair(i,j));
            lastTrace = PTRDIAG;
            ptr -= (Nc+1);
            --i; --j;
         }
      }
      else if (((*ptr)&PTRDIAG) == PTRDIAG && ((*ptr)&PTRLEFT) == PTRLEFT)  { // both left and diag
         if (lastTrace == PTRDIAG) {
            alnidx.push_front(make_pair(i,j));
            ptr -= (Nc+1);
            --i; --j;
         }
         else if (lastTrace == PTRLEFT) {
            alnidx.push_front(make_pair(-1, j));
            --ptr;
            --j;
         }
         else { // choose diag
            alnidx.push_front(make_pair(i,j));
            lastTrace = PTRDIAG;
            ptr -= (Nc+1);
            --i; --j;
         }
      }
      else if (((*ptr)&PTRLEFT) == PTRLEFT && ((*ptr)&PTRTOP) == PTRTOP) { // both Left and Top
         if (lastTrace == PTRLEFT) {
            alnidx.push_front(make_pair(-1, j));
            --ptr;
            --j;
         }
         else if (lastTrace == PTRTOP) {
            alnidx.push_front(make_pair(i, -1));
            ptr -= Nc;
            --i;
         }
         else { // choose left
            alnidx.push_front(make_pair(-1, j));
            lastTrace = PTRLEFT;
            --ptr;
            --j;
         }
      }
      else { // NULL case
         if ((*ptr) != PTRNULL) {
            cerr << "it should be the null case\n";
         }
         break; // for local alnment
      }

      /* old method a little bit simple
      if ( ((*ptr)&PTRDIAG) == PTRDIAG) {
         alnidx.push_front(make_pair(i,j));
         ptr -= (Nc+1);
         --i; --j;
      }
      else if ( ((*ptr)&PTRLEFT) == PTRLEFT) {
         alnidx.push_front(make_pair(-1, j));
         --ptr;
         --j;
      }
      else if ( ((*ptr)&PTRTOP) == PTRTOP) {
         alnidx.push_front(make_pair(i, -1));
         ptr -= Nc;
         --i;
      }
      else if (*ptr == 0) {
         break; // for local alnment
      }
      else {
         debug_showmatrix(cout);
         cerr << "Back pointer: " << *ptr << " at: (i,j) "
            << i << "," << j << " not possible\n";
         exit(1);
      }
      */
   }
   if (alntype == GLOBAL) {
      if (j == -1 && i>0) {
         while (i>-1) {
            alnidx.push_front(make_pair(i--,-1));
         }
      }
      else if (i==-1 && j>0) {
         while (j>-1) {
            alnidx.push_front(make_pair(-1,j--));
         }
      }
   }
}


void Dynaln::findAlnBoundary() {
   if (alntype == LOCAL) {
      seq1begin=alnidx.begin()->first;
      seq2begin=alnidx.begin()->second;
      seq1end=alnidx.rbegin()->first;
      seq2end=alnidx.rbegin()->second;
   }
   else if (alntype == GLOBAL) {
      // For global alignement, the end will be corrected for the terminal gaps.
      seq1begin=0;
      seq2begin=0;
      list<pair<int,int> >::const_iterator it=alnidx.begin();
      while (it != alnidx.end() && it->first == -1) {
         ++seq2begin;
         ++it;
      }
      it=alnidx.begin();
      while (it != alnidx.end() && it->second == -1) {
         ++seq1begin;
         ++it;
      }
      seq1end=seq1->length()-1;
      seq2end=seq2->length()-1;
      list<pair<int,int> >::reverse_iterator li=alnidx.rbegin();
      while (li != alnidx.rend() && li->first == -1) {
         --seq2end;
         ++li;
      }
      li=alnidx.rbegin();
      while (li != alnidx.rend() && li->second == -1) {
         --seq1end;
         ++li;
      }
   }
}

int Dynaln::topAlignBeginIndex() const {
   list<pair<int,int> >::const_iterator it = alnidx.begin();
   while (it != alnidx.end() && it->first == -1) {
      ++it;
   }
   return it->first;
}

int Dynaln::bottomAlignBeginIndex() const {
   list<pair<int,int> >::const_iterator it = alnidx.begin();
   while (it != alnidx.end() && it->second == -1) {
      ++it;
   }
   return it->second;
}

float Dynaln::getNoTerminalGapIdentity() const {
   int e1len = seq1->length() - seq1end - 1;
   int e2len = seq2->length() - seq2end - 1;
   return float(idencnt)/(getAlnlen() - max(seq1begin, seq2begin) - max(e1len, e2len));
}


/* shows the pointer matrix */
void Dynaln::debug_showmatrix(ostream &ous) const {
   int Nc=seq2->length();
   int i,j;
   ous << "\t*";
   for (j=0; j<seq2->length(); j++) ous << (*seq2)[j] << "\t";
   ous << endl;
   for (i=0; i<seq1->length(); i++) {
      ous << (*seq1)[i] << "\t";
      for (j=0; j<seq2->length(); j++) {
         ous << S[i*Nc+j] << "\t";
      }
      ous << endl;
   }
   ous << endl;
}

void Dynaln::showParameters() const {
   cout << "gap insert: " << gapi << " gap extend: " << gape << endl;
   ST.show();
}

// return score
//pair<int,int> Dynaln::global() {
// pointer for traceback (j-i) diag=0, top=1 (Ix), left=-1 (Iy)
// only store the pointers not the scores, too many
int Dynaln::global() throw(AlnInputException) {
   if (seq1->length() < 1) {
      throw AlnInputException("seq1 is empty");
   }
   else if (seq2->length() < 1) {
      throw AlnInputException("seq2 is empty");
   }

   alntype=GLOBAL;
   allocmem();
   int i,j;
   int s1len=seq1->length();  // for short typing
   int s2len=seq2->length();

   M[0]=IX[0]=IY[0]=0;
   M[1]=IX[1]=IY[1]=gapi;
   for (j=2; j<s2len+1; j++) {
      IX[j]=IY[j]=M[j]=M[j-1]+gape;
   }
   int leftM=gapi; 
   int leftIx=gapi;
   int leftIy=gapi;
   int CM, CIx, CIy, backptr; // current values
   int maxDiag, maxTop, maxLeft;
   
#ifdef DEBUG
   cout << "Global alignment matrix:\n";
#endif
   for (i=0; i<s1len; i++) {
#ifdef DEBUG
      cout << "row " << i << endl << leftM << "," << leftIx << "," << leftIy << " | ";
#endif
      for (j=0; j<s2len; j++) {
         maxDiag = max(M[j], max(IX[j], IY[j]));
         CM=maxDiag + ST.lookup(C1[i], C2[j]);
         CIx=max(M[j+1]+gapi, IX[j+1]+gape);
         CIy=max(leftM+gapi, leftIy+gape);
#ifdef DEBUG
         cout << ST.lookup(C1[i], C2[j]) << " " << (*seq1)[i] << "x" << (*seq2)[j] << " "
            << CM << "," << CIx << "," << CIy << " | ";
#endif
         // store pointer
         //Smax=max(CM, max(CIx, CIy));
         if (CM > CIx && CM > CIy) {
            Smax = CM;
            backptr = PTRDIAG;
         }
         else if (CIx > CM && CIx > CIy) {
            Smax = CIx;
            backptr = PTRTOP;
         }
         else if (CIy > CIx && CIy > CM) {
            Smax = CIy;
            backptr = PTRLEFT;
         }
         else if (CM == CIx && CM > CIy) {
            Smax = CM;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            if (maxDiag > maxTop) backptr = PTRDIAG;
            else if (maxTop > maxDiag) backptr = PTRTOP;
            else {
               backptr = (PTRTOP | PTRDIAG);
            }
         }
         else if (CM == CIy && CM > CIx) {
            Smax = CIy;
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxDiag > maxLeft) backptr = PTRDIAG;
            else if (maxLeft > maxDiag) backptr = PTRLEFT;
            else {
               backptr = (PTRLEFT|PTRDIAG); // this is very rare
            }
         }
         else if (CIy == CIx && CIy > CM) {
            Smax = CIy;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxTop > maxLeft) backptr = PTRTOP;
            else if (maxLeft > maxTop) backptr = PTRLEFT;
            else {
               backptr = (PTRLEFT | PTRTOP);
            }
         }
         else { //All three values: CM, CIx, CIy are the same, this is rare!
            Smax = CM;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxDiag > maxTop && maxDiag > maxLeft) backptr= PTRDIAG;
            else if (maxTop > maxDiag && maxTop > maxLeft) backptr= PTRTOP;
            else if (maxLeft > maxDiag && maxLeft > maxTop) backptr= PTRLEFT;
            else if (maxDiag == maxTop && maxDiag > maxLeft) {
               backptr= PTRDIAG | PTRTOP;
            }
            else if (maxDiag == maxLeft && maxDiag > maxTop) {
               backptr= PTRDIAG | PTRLEFT;
            }
            else if (maxTop == maxLeft && maxTop > maxDiag) {
               backptr= PTRTOP | PTRLEFT;
            }
            else {
               cerr << "it is very rare that the three values are identical: "
                  << maxDiag << ", " << maxTop << ", " << maxLeft << endl;
               backptr= PTRDIAG | PTRTOP | PTRLEFT;
            }
         }
         //cout << backptr << "\t";
         S[i*s2len+j] = backptr;
         // ready for the next round
         M[j]=leftM; IX[j]=leftIx; IY[j]=leftIy;
         leftM=CM; leftIx=CIx; leftIy=CIy;
      }
#ifdef DEBUG
      cout << endl;
#endif
      M[j]=CM; IX[j]=CIx; IY[j]=CIy;
      leftM=leftIx=leftIy=(M[0] + gape);
   }
   Smaxi=seq1->length()-1;
   Smaxj=seq2->length()-1;
   clearResult();
   //cerr << "global score: " << Smax << endl;
   return Smax;
}

int Dynaln::local() {
   alntype=LOCAL;
   allocmem();
   int i,j,s1len, s2len;
   s2len=seq2->length();
   s1len=seq1->length();
   Smax=Smaxi=Smaxj=0;
#ifdef DEBUG
   cerr << "\ngap parameter inside local(): " << gapi << " " << gape << endl;
#endif

   // fill the max array with zero
   // we only need the top row for memory
   // 0 as bourder value, 1 as first value, seq2len as last value
   for (j=0; j <= s2len; j++) M[j]=IX[j]=IY[j]=0;
   int leftM=0; 
   int leftIx=0;
   int leftIy=0;
   int CM, CIx, CIy, backptr, currScore; // current values
   /*
    * Use 3 arrays, plus left M,Ix,Iy to represent the 3 matrices: match, insert gap in X,
    * insert gap in Y.
    * X for column , Y row
    * o---- Y ---->
    * |
    * X
    * |
    *\|/
    * Here I have a little bit of code duplication for performance.
    * Use CM, CIx, CIy to decide the trace back pointer direction.
    * If there is a tie, then use maxima of diag, top, and left
    * to decide which direction to point to. Even with this
    * two layers of inforamtion, there is still cases where the
    * scores are identical for all three of them. The pointer is
    * saved as an OR operation, the it is left to the trace-back
    * algorithm to decide which direction to pick.
    */
   
   int maxDiag, maxTop, maxLeft;
   for (i=0; i < s1len; i++) {
#ifdef DEBUG
      cout << "row i " << i << " " << (*seq1)[i] << endl;
#endif
      for (j=0; j<s2len; j++) {
         //cout << "col j " << j << " " << (*seq2)[j];
#ifdef DEBUG
         cout << " " << j << (char)toupper((*seq2)[j]);
#endif
         // j is the top row previus cell, C2[j] is the current sequence
         maxDiag = max(M[j], max(IX[j], IY[j]));
         CM=maxDiag + ST.lookup(C1[i], C2[j]);
         CIx=max(M[j+1]+gapi, IX[j+1]+gape); // above cell
         CIy=max(leftM+gapi, leftIy+gape);
         if (CM <= 0 && CIx <= 0 && CIy <= 0) {
            currScore = 0;
            backptr = PTRNULL;
         }
         else if (CM > CIx && CM > CIy) {
            currScore = CM;
            backptr = PTRDIAG;
         }
         else if (CIx > CM && CIx > CIy) {
            currScore = CIx;
            backptr = PTRTOP; // insert gap in seq2
         }
         else if (CIy > CM && CIy > CIx) {
            currScore = CIy;
            backptr = PTRLEFT; // insert gap in seq1
         }
         else if (CM == CIx && CM > CIy) {
            currScore = CIx;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            if (maxDiag > maxTop) backptr = PTRDIAG;
            else if (maxTop > maxDiag) backptr = PTRTOP;
            else {
               //cerr << i << "," << j 
               //   << " Top and Diag trace pointer have the same probability!\n";
               backptr = (PTRTOP | PTRDIAG);
            }
         }
         else if (CM == CIy && CM > CIx) {
            currScore = CIy;
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxDiag > maxLeft) backptr = PTRDIAG;
            else if (maxLeft > maxDiag) backptr = PTRLEFT;
            else {
               //cerr << i << ", " << j 
               //   << " Left and Diag trace pointer have the same probability!\n";
               backptr = (PTRLEFT|PTRDIAG); // this is very rare
            }
         }
         else if (CIy == CIx && CIy > CM) {
            currScore = CIy;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            if (maxTop > maxLeft) backptr = PTRTOP;
            else if (maxLeft > maxTop) backptr = PTRLEFT;
            else {
               //cerr << i << ", " << j
               //   << " Left and Top track pointer have the same probability!\n";
               backptr = (PTRLEFT | PTRTOP);
            }
         }
         else { // this state should be impossible!
            currScore = CM;
            maxTop = max(M[j+1], max(IX[j+1], IY[j+1]));
            maxLeft = max(leftM, max(leftIy, leftIx));
            //cerr << i << ", " << j
            //   << " All three trace pointer have the same probability! very rare\n"
            //   << " Scores for CM, CIx, CIy: "
            //   << CM << ", " << CIx << ", " << CIy << endl
            //   << " Diag, Top, Left: " << maxDiag << ", " << maxTop << ", "
            //   << maxLeft << endl;
            if (maxDiag > maxTop && maxDiag > maxLeft) backptr= PTRDIAG;
            else if (maxTop > maxDiag && maxTop > maxLeft) backptr= PTRTOP;
            else if (maxLeft > maxDiag && maxLeft > maxTop) backptr= PTRLEFT;
            else if (maxDiag == maxTop && maxDiag > maxLeft) {
               //cerr << "Diag and Top trace pointer have the same probability!\n";
               backptr= PTRDIAG | PTRTOP;
            }
            else if (maxDiag == maxLeft && maxDiag > maxTop) {
               //cerr << "Diag and Left trace pointer have the same probability!\n";
               backptr= PTRDIAG | PTRLEFT;
            }
            else if (maxTop == maxLeft && maxTop > maxDiag) {
               //cerr << "Top and Left trace pointer have the same probability!\n";
               backptr= PTRTOP | PTRLEFT;
            }
            else {
               cerr << "it is very rare that the three values are identical: "
                  << maxDiag << ", " << maxTop << ", " << maxLeft << endl;
               backptr= PTRDIAG | PTRTOP | PTRLEFT;
            }
         }
#ifdef DEBUG
         cout << " (" << CM << ", " << CIx << ", " << CIy << ") ";
#endif
         // save the trace-back pointer in the matrix
         S[i*s2len+j] = backptr;
         if (currScore > Smax) {
            Smax=currScore;
            Smaxi=i; Smaxj=j;
         }
         // ready for the next round
         M[j]=leftM; IX[j]=leftIx; IY[j]=leftIy; // for the lower row
         leftM=CM; leftIx=CIx; leftIy=CIy;
      }
      //cout << " max: " << Smax << " [" << Smaxi << ", " << Smaxj << "]\n";
      M[j]=CM; IX[j]=CIx; IY[j]=CIy;
      leftM=leftIx=leftIy=0;
   }
   clearResult();
#ifdef DEBUG
   // show pointer matrix
   debug_showmatrix(cout);
   cout << "max score " << Smax << endl;
#endif
   return Smax;
}


void Dynaln::countnumgaps() {
   list<pair<int,int> >::const_iterator li=alnidx.begin();
   numgaps1=0;
   numgaps2=0;
   while (li != alnidx.end()) {
      if (li->first == -1) {
         ++numgaps1;
         while (li != alnidx.end() && li->first == -1) ++li;
      }
      else ++li;
   }
   li=alnidx.begin();
   while (li != alnidx.end()) {
      if (li->second == -1) {
         ++numgaps2;
         while (li != alnidx.end() && li->second == -1) ++li;
      }
      else ++li;
   }
   //cerr << "gaps: " << numgaps1 << " " << numgaps2 << endl;
}

void Dynaln::buildResult(const int delta1, const int delta2) {
   //cerr << "trace pointer ...\n";
   //traceback(Smaxi, Smaxj);
   tracepointer(Smaxi, Smaxj);
   //clearResult(); // this is called just after calling global() or local()
   //cerr << "find aln boundary ...\n";
   findAlnBoundary();
   /* reset by clearResult()
   idencnt=0;
   simcnt=0;
   numgaps1=0;
   numgaps2=0;
   topaln.clear();
   bottomaln.clear();
   middle.clear();
   */
   //cerr << "build aln info ...\n";
   buildAlnInfo(delta1, delta2);
   //cerr << "count number of gaps ...\n";
   countnumgaps();
   //cerr << "buildResult done\n";
}

/*
 * marking at every 10 this could become a parameter if you
 * use very long sequence to do the alignment
 */
void Dynaln::buildAlnInfo(const int delta1, const int delta2) {
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i,j,counter;
   counter=0;
   char topchar, bottomchar;
   vector<int> topmark(alnidx.size(),-1);
   vector<int> bottommark(alnidx.size(),-1);

   while (lit != alnidx.end()) {
      i=lit->first; // top sequence index
      j=lit->second; // bottom sequence index
      //if (counter % 10 == 0) {
      if (counter % markevery == 0) {
         if (i>-1) { 
            topmark[counter]=i+1+delta1;
         }
         if (j>-1) {
            bottommark[counter]=j+1+delta2;
         }
      }
      if (i>=0) {
         topchar=toupper((*seq1)[i]); 
         //topchar=(*seq1)[i];
      }
      else {
         topchar=gapchar;
         ++gaplen1;
      }
      topaln += topchar;

      if (j>=0) {
         //bottomchar = (*seq2)[j];
         bottomchar = toupper((*seq2)[j]);
      }
      else {
         bottomchar = gapchar;
         ++gaplen2;
      }
      bottomaln += bottomchar;

      if (topchar != gapchar && bottomchar != gapchar) {
         if (topchar == bottomchar) {
            middle += idenchar;
            ++idencnt;
         }
         else if (ST.lookup(C1[i],C2[j]) >= simcut) {
            ++simcnt;
            middle += simchar;
         }
         else middle += ' ';
      }
      else {
         middle += ' ';
      }
      ++lit; ++counter;
   }
   string topr(alnidx.size(), ' ');
   string bottomr(alnidx.size(), ' ');
   string x;

   for (i=0; i<alnidx.size(); i++) {
      if (topmark[i] > -1) {
         x=itos(topmark[i]);
         topr.replace(i, x.length(), x);
      }
      if (bottommark[i] > -1) {
         x=itos(bottommark[i]);
         bottomr.replace(i, x.length(), x);
      }
   }
   topruler=topr;
   bottomruler=bottomr;
}

// the implementations are similar to getNgMatchArray
// I did not make a one depends on the other. The performance
// will be slower that way.
void Dynaln::getMatchArray(vector<int> &marr, bool countsim) const {
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i,j,k,count;
   //vector<int> marr(alnidx.size(), 0);
   marr.resize(alnidx.size());
   k=0;
   while (lit != alnidx.end()) {
      i=lit->first;
      j=lit->second;
      if (i<0 || j<0) count= - 1;
      else if ((*seq1)[i] == (*seq2)[j]) count=2;
      else if (countsim && ST.lookup(C1[i],C2[j]) >= simcut) count=1;
      else count=0;
      marr[k++]=count;
      ++lit;
   }
}

// this version only cares about the non-gapped version
// this is being used for bootstrap
// this is used by bootstrap
void Dynaln::getNgMatchArray(vector<int> &ngmarr, bool append, bool countsim ) const {
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i, j, count;
   if (!append) ngmarr.clear();
   while (lit != alnidx.end()) {
      i=lit->first;
      j=lit->second;
      if (i>=0 && j>=0) {
         if ((*seq1)[i] == (*seq2)[j]) count=2;
         else if (countsim && ST.lookup(C1[i],C2[j]) >= simcut) count=1;
         else count=0;
         ngmarr.push_back(count);
      }
      ++lit;
   }
}

void Dynaln::printAlign(ostream &ous, const int w) const {
   //cerr << "print align info for human ...\n";
   ous << seq1->getName() << " x " << seq2->getName() << "  "
      << seq1begin << "-" << seq1end << "/" << seq1->length()
      << " | " 
      << seq2begin << "-" << seq2end << "/" << seq2->length()
      << endl
      << "Score=" << Smax 
      << " gap length: " << gaplen1 << " " << gaplen2
      << " num gaps: " << numgaps1 << " " << numgaps2
      << " idencnt=" << idencnt 
      << " simcnt=" << simcnt
      << " alnlen=" << topaln.length()
      << " identity=" << setprecision(5) << getIdentity()
      << " similarity=" << getSimilarity()
      << endl;
   // topaln, middle and bottomaln should be the same length
   if (topaln.length() != middle.length() 
         || middle.length() !=bottomaln.length()) {
      cerr << "the final alignment string is not right\n";
      exit(1);
   }
   //int i=0, i1=1, i2=1, charcnt;
   int i=0, j;
   //string tmp, tmp2, ruler1, ruler2;
   string ruler(w,' '), ruler1, ruler2;
   for (int x=0; x<w; x+=markevery) { // markevery default 10 residues
      ruler[x]='+';
   }
   while (i<topaln.length()) {
      j=i+w-1;
      while (isdigit(topruler[j])) ++j;
      ruler1=topruler.substr(i, j-i);
      j=i+w-1;
      while (isdigit(bottomruler[j])) ++j;
      ruler2=bottomruler.substr(i, j-i);

      ous << ruler1 << endl << ruler << endl
         << topaln.substr(i,w) << endl
         << middle.substr(i, w) << endl
         << bottomaln.substr(i,w) << endl << ruler << endl
         << ruler2 << endl << endl;
      i += w;
   }
}


pair<string,int> Dynaln::getNucleicConsensus() const {
   string tmp(bottomaln);
   int diff = 0;
   for (size_t i = 0; i<tmp.size(); ++i) {
      if (tmp[i] != topaln[i]) {
         if (tmp[i] == gapchar) tmp[i] = topaln[i];
         else {
            if ((tmp[i] == 'A' && topaln[i] == 'C') ||
                  (tmp[i] == 'C' && topaln[i] == 'A')) tmp[i] = 'M';
            else if ((tmp[i] == 'A' && topaln[i] == 'G')
                  || (tmp[i] == 'G' && topaln[i] == 'A')) tmp[i] = 'R';
            else if ((tmp[i] == 'A' && topaln[i] == 'T') ||
                  (tmp[i] == 'T' && topaln[i] == 'A')) tmp[i] = 'W';
            else if ((tmp[i] == 'C' && topaln[i] == 'G') ||
                  (tmp[i] == 'G' && topaln[i] == 'C')) tmp[i] = 'S';
            else if ((tmp[i] == 'C' && topaln[i] == 'T') ||
                  (tmp[i] == 'T' && topaln[i] == 'C')) tmp[i] = 'Y';
            else if ((tmp[i] == 'G' && topaln[i] == 'T') ||
                  (tmp[i] == 'T' && topaln[i] == 'G')) tmp[i] = 'K';
            ++diff;
         }
      }
   }
   return pair<string, int>(tmp, diff);
}

pair<string,string> Dynaln::getNucleicConsensus2() const {
   string tmpt(topaln);
   string tmpb(bottomaln);
   int diff=0;
   for (size_t i = 0; i<tmpt.size(); ++i) {
      if (tmpt[i] == gapchar) tmpt[i] = tmpb[i];
      else if (tmpb[i] == gapchar) tmpb[i] = tmpt[i];
      else if (tmpt[i] != tmpb[i]) {
         ++diff;
      }
   }
   if (diff == 0) {
      return make_pair(tmpt, string());
   }
   else {
      return pair<string, string>(tmpt, tmpb);
   }
}

void Dynaln::clearResult() {
   if (idencnt == 0) return;
   idencnt=0;
   simcnt=0;
   numgaps1=0;
   numgaps2=0;
   gaplen1=0;
   gaplen2=0;
   //alnlen=0;
   topaln.clear();
   middle.clear();
   bottomaln.clear();
}

string Dynaln::toDelimitedString(const string &dl, int ibase) const {
   ostringstream ous;
   try {
      /* production version
      ous << seq1->getName() << dl << seq1->length() << dl
         << seq2->getName() << dl << seq2->length() << dl
         << Smax << dl << idencnt << dl << simcnt << dl
         << alnidx.size() << dl
         << numgaps1 << dl << numgaps2 << dl
         << gaplen1 << dl << gaplen2 << dl
         << seq1begin << dl <<  seq1end << dl
         << seq2begin << dl <<  seq2end << dl
         << setprecision(4)
         << (seq1->computeEntropy(seq1begin, seq1end)).first << dl
         << (seq2->computeEntropy(seq2begin, seq2end)).first << dl;
      */
      // debug version, so we know which line crashed
      //pair<double,double> tmp=seq1->computeEntropy();
      //cerr << tmp.first << " " << tmp.second << endl;
      //cout << seq1begin << " " << seq1end << endl;
      ous << seq1->getName() << dl << seq1->length() << dl;
      ous   << seq2->getName() << dl << seq2->length() << dl;
      ous   << Smax << dl << idencnt << dl << simcnt << dl;
      ous   << alnidx.size() << dl;
      ous   << numgaps1 << dl << numgaps2 << dl;
      ous   << gaplen1 << dl << gaplen2 << dl;
      ous   << seq1begin+ibase << dl <<  seq1end+ibase << dl;
      ous   << seq2begin+ibase << dl <<  seq2end+ibase << dl;
      ous << setprecision(4) << (seq1->computeEntropy(seq1begin, seq1end)).first << dl;
      ous   << (seq2->computeEntropy(seq2begin, seq2end)).first << dl;
   }
   catch (exception &er) {
      cerr << er.what() << endl;
      cerr << "fetal error inside Dynaln::toDelimitedString()" << endl;
      exit(1);
   }
   return ous.str();
}

string Dynaln::headers(const string &dl) {
   return "seq1name" + dl + "seq1len" + dl
      + "seq2name" + dl + "seq2len" + dl
      + "score" + dl
      + "identical" + dl + "similar" + dl + "alnlen" + dl
      + "seq1numgap" + dl + "seq2numgap" + dl
      + "seq1gaplen" + dl + "seq2gaplen" + dl
      + "seq1begin" + dl + "seq1end" + dl
      + "seq2begin" + dl + "seq2end" + dl
      + "seq1entropy" + dl + "seq2entropy";
}

int Dynaln::getSeq1Length() const throw (runtime_error) {
   if (seq1 == 0) throw runtime_error("NULL pointer for aligner's seq1");
   return seq1->length();
}

int Dynaln::getSeq2Length() const throw (runtime_error) {
   if (seq2 == 0) throw runtime_error("NULL pointer for aligner's seq2");
   return seq2->length();
}

////////////////// Linear Space Algorithms ///////////////////////////////

// return the pointers
//void Dynaln::allocmemLS() {
void LSDynaln::allocmem() {
   int Nr=seq1->length()+1;
   int Nc=seq2->length()+1;
   if (M == 0) {
      Ssize=2*Nc;
      numcol=Nc;
      //S = new int[Nc];
      M=new int[Ssize];
      IX = new int[Ssize];
      IY = new int[Ssize];
      //Itop = new int[Nc];  // top insert Y array
      //SR = new int[Nc];
      MR=new int[Ssize];
      IXR = new int[Ssize];
      IYR = new int[Ssize];
      //if (alntype == GLOBAL) {
        // ItopR = new int[Nc];
      //}
   }
   else {
      if (Ssize < 2*Nc) {
         //delete[] S; 
         delete[] M; delete[] IX; delete[] IY;
         //delete[] SR; 
         delete[] MR; delete[] IXR; delete[] IYR;

         Ssize=2*Nc;
         numcol=Nc;

         //S = new int[Nc];
         M=new int[Ssize];
         IX = new int[Ssize];
         IY = new int[Ssize];
         //Itop = new int[Nc];  // top insert Y array
         //SR = new int[Nc];
         MR=new int[Ssize];
         IXR = new int[Ssize];
         IYR = new int[Ssize];
      }
   }
   C1=seq1->getcode();
   C2=seq2->getcode();
}

// helper debug function
void showMatrix(pair<int*,int*> &M, int w) {
   for (int i=0; i<w; i++) {
      cerr << M.first[i] << " ";
   }
   cerr << endl;
   for (int i=0; i<w; i++) {
      cerr << M.second[i] << " ";
   }
   cerr << endl;
}
      
/* the correct version, using 4 arrays
 * M, IX, IY, S, and 4 for reverse
 * Separating the two algorithms, backward and forward
 * return the start index of top and bottom rows 
 * respectively.
 */
pair<int,int> LSDynaln::computeScoreFW(int b1, int e1, int b2, int e2) 
{
   int Nc=e2-b2+2; // number of the clumns of the scoring matrix
   //int Nr=e1-b1+2;
   int i,j;
   M[0]=IX[0]=IY[0]=0;
   M[1]=IX[1]=IY[1]=gapi;
   for (j=2; j<Nc; j++) {
      IY[j]=IX[j]=M[j]=M[j-1]+gape;
   }
   M[Nc]=IX[Nc]=IY[Nc]=gapi;
   // matrix is 1 more than the sequence length
   int *MDia, *MLeft, *IXTop, *IYDia, *IXLeft, *IYLeft;
   int toprow=0;
   int bottomrow=Nc;
 /*
#ifdef DEBUG
   cout << "Matrix forward from computeScoreFW " << b1 << "-" << e1
      << " x " << b2 << "-" << e2 << "\n";
#endif
*/
   for (i=b1; i<=e1; i++) {
      MDia=M+toprow;  // MTop = MDia+1, no need to keep
      IYDia=IY+toprow; 
      IXTop=IX+toprow+1; //IXDia=IX+toprow; not needs
      MLeft=M+bottomrow;
      IXLeft=IX+bottomrow;
      IYLeft=IY+bottomrow;
      /*
#ifdef DEBUG
      cout << "row " << i-b1 << endl << *MLeft << "," << *IXLeft << "," << *IYLeft << " | ";
#endif
*/
      for (j=b2; j<=e2; j++) {
         *(MLeft+1)=max(*MDia, max(*(IXTop-1), *IYDia)) + ST.lookup(C1[i], C2[j]);
         *(IXLeft+1)=max(*(MDia+1)+gapi, *IXTop+gape);
         *(IYLeft+1)=max(*MLeft+gapi, *IYLeft+gape);
         /*
#ifdef DEBUG
         cout << ST.lookup(C1[i],C2[j]) << " " << (*seq1)[i] << "x" << (*seq2)[j] << " "
            << *(MLeft+1) << "," << *(IXLeft+1) << "," << *(IYLeft+1) << " | ";
#endif
*/
         ++MDia; ++IXTop; ++IYDia;
         ++MLeft; ++IXLeft; ++IYLeft;
      }
      /*
#ifdef DEBUG
      cout << endl;
#endif
*/
      // initialize the first column, boundary condition
      // This overwrites the top row first column if done
      M[toprow]=IX[toprow]=IY[toprow]=M[bottomrow]+gape;
      toprow=Nc-toprow;
      bottomrow=Nc-bottomrow;
   }
   // reverse the first column initialization effect
   // now the top and bottom has also been switched
   M[bottomrow]=IX[bottomrow]=IY[bottomrow]=M[toprow]-gape;
   //cout << "====================================\n";
   //return make_pair(toprow, bottomrow);
   //need to reverse the last operation
   return make_pair(bottomrow, toprow);
}
// (b1>e1 && b2>=e2) || (b1>=e1 && b2>e2)
pair<int,int> LSDynaln::computeScoreBW(int b1, int e1, int b2, int e2) 
{
   int Nc=b2-e2+2;
   //int Nr=b1-e1+2;
   int i,j;
   MR[0]=IXR[0]=IYR[0]=0;
   MR[1]=IXR[1]=IYR[1]=gapi;
   for (j=2; j<Nc; j++) {
      IYR[j]=IXR[j]=MR[j]=MR[j-1]+gape;
   }
   MR[Nc]=IXR[Nc]=IYR[Nc]=gapi;
   // matrix is 1 more than the sequence length
   int *MDia, *MLeft, *IXTop, *IYDia, *IXLeft, *IYLeft;
   //int toprow, bottomrow, r;
   int toprow=0;
   int bottomrow=Nc;
   /*
#ifdef DEBUG
   cout << "Matrix bacdward from computeScoreBW " << b1 << "-" << e1
      << " x " << b2 << "-" << e2 << "\n";
   for (j=0; j<Nc; j++) {
      cout << MR[j] << "," << IXR[j] << "," << IYR[j] << " | ";
   }
#endif
*/
   for (i=b1; i>=e1; i--) {
      MDia=MR+toprow;  // MTop = MDia+1, no need to keep
      IXTop=IXR+toprow+1; //IXDia=IX+toprow; not needs
      IYDia=IYR+toprow;
      MLeft=MR+bottomrow;
      IXLeft=IXR+bottomrow;
      IYLeft=IYR+bottomrow;
      //cout << "row " << b1-i << endl << *MLeft << ", " << *IXLeft << ", " << *IYLeft << " | ";
      for (j=b2; j>=e2; j--) {
         *(MLeft+1)=max(*MDia, max(*(IXTop-1), *IYDia)) + ST.lookup(C1[i], C2[j]);
         *(IXLeft+1)=max(*(MDia+1)+gapi, *IXTop+gape);
         *(IYLeft+1)=max(*MLeft+gapi, *IYLeft+gape);
         //copying and  ready for the next round
         /*
#ifdef DEBUG
         cout << ST.lookup(C1[i], C2[j]) << " " << (*seq1)[i] << "x" << (*seq2)[j] << " "
            << *(MLeft+1) << "," << *(IXLeft+1) << "," << *(IYLeft+1) << " | ";
#endif
*/
         ++MDia; ++IXTop; ++IYDia;
         ++MLeft; ++IXLeft; ++IYLeft;
      }
      //cout << endl;
      // initialize the first column, boundary condition
      MR[toprow]=IXR[toprow]=IYR[toprow]=MR[bottomrow]+gape;
      toprow=Nc-toprow;
      bottomrow=Nc-bottomrow;
   }
   // reverse the intialization effect (used for looping setup) when done
   MR[bottomrow]=IXR[bottomrow]=IYR[bottomrow]=MR[toprow]-gape;
   return make_pair(bottomrow, toprow);
}

int LSDynaln::global() throw(AlnInputException) {
   if (seq1->length() < 1) {
      throw AlnInputException("seq1 empty");
   }
   else if (seq2->length() < 1) {
      throw AlnInputException("seq1 empty");
   }

   cerr << "\n*** Global linear space dynamic alignment ***\n";
   alntype=GLOBAL;
   allocmem();
   if (!alnidx.empty()) {
      alnidx.clear();
   }
   Smax=path(0,seq1->length()-1, 0, seq2->length()-1);
   clearResult();
}

// b1 < e1 && b2 < e2
int LSDynaln::path(int b1, int e1, int b2, int e2) {
   int row=e1-b1+1; // length of seq1 range 
   int col=e2-b2+1;  // length of seq range
   pair<int, int> forward, backward;
   list<pair<int,int> > pathForward, pathBackward;
   list<pair<int,int> >::const_iterator it;
   int s1begin, s1end, s2begin, s2end; // divide boundary
   int scoreFB;
   //int i1,j1, i2,j2;
#ifdef DEBUG
   cout << "\nfinding best path for  " << b1 << "-" << e1
      << " x " << b2 << "-" << e2 << endl;
      //<< seq1->substr(b1+1, e1+1) << endl
      //<< seq2->substr(b2+1, e2+1) << endl;
#endif

   if (row < 2 && col < 2) {
      //cout << b1 << "," << b2 << " |\n";
      alnidx.push_back(make_pair(b1,b2));
      return ST.lookup((*seq1)[b1], (*seq2)[b2]);
   }
   else if (row < 2 || col < 2) {
      forward=computeScoreFW(b1,e1,b2,e2);
      s1end=e1;
      s2end=e2;
      pathForward=tracebackFW(forward,b1,s1end,b2,s2end);
      alnidx.insert(alnidx.end(), pathForward.begin(), pathForward.end());
      int bb=forward.second;
      return max(M[bb+col], max(IX[bb+col], IY[bb+col]));
   }

   s1end=(b1+e1)/2;  // top part end
   s1begin=s1end+1; // bottom part begin
   forward=computeScoreFW(b1, s1end, b2, e2);
   backward=computeScoreBW(e1, s1begin, e2, b2);
   //int *F=S+forward.second;
   //int *B=S+backward.second;
   ///////// forward ///////////
   int *scoreM=M+forward.second;
   int *scoreIX=IX+forward.second;
   int *scoreIY=IY+forward.second;
   ////// backward //////////////
   int *scoreMR=MR+backward.second;
   int *scoreIXR=IXR+backward.second;
   int *scoreIYR=IYR+backward.second;

   int maxScore=-99999999; // a very small number
   int maxj, j;
   maxj=-1;
   //cout << "find the division cell\n";
   for (j=0; j<=col; j++) {
      scoreFB=max(scoreM[j], max(scoreIX[j], scoreIY[j]))
            + max(scoreMR[col-j], max(scoreIXR[col-j], scoreIYR[col-j]));
      /*
#ifdef DEBUG
      cout << scoreFB << " " << max(scoreM[j], max(scoreIX[j], scoreIY[j]))
            << " + " << max(scoreMR[col-j], max(scoreIXR[col-j], scoreIYR[col-j]))
            << " | ";
#endif
*/
      //if (F[j] + B[col-j] > maxScore) {
      if (scoreFB > maxScore) {
         maxScore= scoreFB;
         maxj=j;
      }
   }
   //cout << endl;
   // divide the large problem into the following two smaller ones
   // The index is in terms of the sequences
   // seq1 b1----s1end  s1begin----e1
   //      |||||||||||  |||||||||||||
   // seq2 b2----s2end  s2begin----e2
   // The traceback functions works on indices on the matrix which
   // has an additional dimention
   /*
   if (maxj == 0) {
      cerr << "Abnormal division result with maxscore: " << maxScore
         << " and maxj: " << maxj << endl
         << seq1->substr(b1+1, e1+1) << "\n" << seq2->substr(b2+1,e2+1)
         << endl << "seq1 divisions " << s1end << " | " << s1begin << endl;
      debug_showmatrixForward(forward,col+1);
      debug_showmatrixBackward(backward,col+1);
      exit(1);
   }
   */
   s2end=b2+maxj-1;
   s2begin=e2-col+maxj+1;
#ifdef DEBUG
   cout << "new seq1 division point: " << s1end << " " << s1begin << endl;
   cout << "new seq2 division point: " << s2end << " " << s2begin << endl;
   // speed up by skipping identical residues
#endif
   bool topDone=false;
   if (maxj == 0)  {
      topDone=true;
      while (s1end >= b1) {
         pathForward.push_front(make_pair(s1end--, -1));
      }
   }
   else if (s2end == b2 && b1 == s1end) {
      topDone=true;
      pathForward.push_front(make_pair(b1,b2));
   }
   else {
      //cout << "Before tracebackFW s1end,s2end(" << s1end << ", " << s2end << ")\n";
      pathForward = tracebackFW(forward, b1, s1end, b2, s2end);
      //cout << "After tracebackFW s1end,s2end(" << s1end << ", " << s2end << ")\n";
      if (s1end <= b1 || s2end <= b2) {
         //cout << "the new ends indicates finish of job after tracebackFW! "
            //<< b1 << "--" << s1end << " |x| " 
            //<< b2 << "--" << s2end << "\n";
         topDone=true;
         //exit(1);
      }
   }
#ifdef DEBUG
   cerr << "the ends before tracebackBW! "
      << s1begin << "--" << e1 << " || " 
      << s2begin << "--" << e2 << "\n";
   if (s2begin > e2) {
      cerr << "Division finished the bottom problem with maxscore: " << maxScore
         << " and maxj: " << maxj << endl
         << seq1->substr(b1+1, e1+1) << "\n" << seq2->substr(b2+1,e2+1)
         << endl << "seq1 divisions " << s1end << " | " << s1begin << endl;
      debug_showmatrixForward(forward,col+1);
      debug_showmatrixBackward(backward,col+1);
   }
#endif
   bool bottomDone=false;
   if (maxj == col) { // s2begin == e2+1
      bottomDone=true;
      while (s1begin <= e1) {
         pathBackward.push_back(make_pair(s1begin++, -1));
      }
   }
   else if (s2begin == e2 && e1==s1begin) {
      bottomDone=true;
      pathBackward.push_back(make_pair(e1,e2));
   }
   else {
      pathBackward = tracebackBW(backward, e1, s1begin, e2, s2begin);
      if (s1begin >= e1 || s2begin >= e2) {
         //cerr << "new ends indicated finish of the bottom job after tracebackBW! "
            //<< s1begin << "--" << e1 << " || " 
            //<< s2begin << "--" << e2 << "\n";
         bottomDone=true;
         //exit(1);
      }
   }
   pathForward.insert(pathForward.end(), pathBackward.begin(), pathBackward.end());
   if (!topDone) path(b1, s1end, b2, s2end);
   // save the result into the shared alnidx list
   alnidx.insert(alnidx.end(), pathForward.begin(), pathForward.end());
   if (!bottomDone) path(s1begin, e1, s2begin, e2);
   return maxScore;
}

/* forward version, b1<e1 && b2<e2
 * the matrix has one additional row,column than the sequences
 */
list<pair<int,int> > LSDynaln::tracebackFW(pair<int, int> mrow, int b1, int &e1, int b2, int &e2) 
{
   list<pair<int, int> > path;
   int *MTop=M+mrow.first;
   int *MButtom=M+mrow.second; 
   int *IXTop=IX+mrow.first; 
   int *IXButtom=IX+mrow.second;
   int *IYTop=IY+mrow.first;
   int *IYButtom=IY+mrow.second;
   int r,c;
   r=e1-b1+1;
   c=e2-b2+1; // column index in the matrix

   int R=r;  // remember the last row
   //int i,j; // index in sequence for reverse lookup
   // stop the loop after moving into the top row
   //i=e1;
   //j=e2;
   int score;
   // stop when moved to the top row.
   while (c>0 && R-r != 1 && e2 >= b2) {
      score=max(MButtom[c], max(IXButtom[c], IYButtom[c]));
      if (score == MButtom[c]) {
         path.push_front(make_pair(e1, e2));
         --r; --c;
         --e1; --e2;
      }
      else if (score == IXButtom[c]) {
         path.push_front(make_pair(e1, -1));
         --r; --e1;
      }
      else if (score == IYButtom[c]) { 
         // from left, trace till it ends
         path.push_front(make_pair(-1, e2));
         --c; --e2;
      }
      else {
         cerr << "track back error. not possible\n";
         exit(1);
      }
   }
   // we are not working in the top row, and only move to the left
   while (c>0 && R-r == 1 && e2>=b2) {
      score=max(MTop[c], max(IXTop[c], IYTop[c]));
      if (score==MTop[c]) {
         path.push_front(make_pair(e1, e2));
         --r; --c; --e1; --e2;
      }
      else if (score == IXTop[c]) {
         path.push_front(make_pair(e1, -1));
         --r; --e1;
      }
      else if (score == IYTop[c]) {
         path.push_front(make_pair(-1, e2));
         --c; --e2;
      }
   }
   while (c==0 && r>0 && e1 >= b1) {
      path.push_front(make_pair(e1--, -1));
      --r;
   }
   //if (i>=b1) e1=i;
   //if (j>=b2) e2=j;
   return path;
}

// backward version, b1 >= e1 && b2 >= e2
// at this point having two separate version is better than keeping one for both
// Essentially copying and pasting
list<pair<int,int> > LSDynaln::tracebackBW(pair<int, int> mrow, int b1, int &e1, int b2, int &e2) 
{
   list<pair<int, int> > path;
   int *MTop=M+mrow.first;
   int *MButtom=M+mrow.second; 
   int *IXTop=IX+mrow.first; 
   int *IXButtom=IX+mrow.second;
   int *IYTop=IY+mrow.first;
   int *IYButtom=IY+mrow.second;
   int r,c;
   r=b1-e1+1;
   c=b2-e2+1;

   int R=r;  // remember the last row
   //int i,j;  // index in sequence for reverse lookup
   // stop the loop after moving into the top row
   //i=e1;
   //j=e2;
   int score;
   // stop when moved to the top row.
   while (c>0 && R-r < 1 && e2 <= b2) {
      score=max(MButtom[c], max(IXButtom[c], IYButtom[c]));
      if (score == MButtom[c]) {
         path.push_back(make_pair(e1, e2));
         --r; --c; ++e1; ++e2;
      }
      else if (score == IXButtom[c]) {
         path.push_back(make_pair(e1, -1));
         --r; ++e1;
      }
      else if (score == IYButtom[c]) { 
         // from left, trace till it ends
         path.push_back(make_pair(-1, e2));
         --c; ++e2;
      }
      else {
         cerr << "track back error. not possible\n";
         exit(1);
      }
   }
   // we are not working in the top row, and only move to the left
   while (c>0 && R-r == 1 && e2 <= b2) {
      score=max(MTop[c], max(IXTop[c], IYTop[c]));
      if (score==MTop[c]) {
         path.push_back(make_pair(e1, e2));
         --r; --c; ++e1; ++e2;
      }
      else if (score == IXTop[c]) {
         path.push_back(make_pair(e1, -1));
         --r; ++e1;
      }
      else if (score == IYTop[c]) {
         path.push_back(make_pair(-1, e2));
         --c; ++e2;
      }
   }
   // special case trace when column 0 of the matrix is reached
   while (c==0 && r>0 && e1 <= b1) {
      path.push_back(make_pair(e1++, -1));
      --r;
   }
   //if (i<=b1) e1=i;
   //if (j<=b2) e2=j;
   return path;
}

// overwrite the square version
void LSDynaln::buildResult(const int delta1, const int delta2) {
   //cerr << " *** calling LSDynaln buildResult() \n";
   list<pair<int, int> >::const_iterator lit= alnidx.begin();
   int i,j;
   char topchar, bottomchar;
   //clearResult();
   // find sequence begin and end, removing intial gaps
   // The local version produce this result in the
   // matrix computation step.
   if (alntype==GLOBAL) {
      seq1begin=seq2begin=0;
      seq1end=seq1->length()-1;
      seq2end=seq2->length()-1;
      if (lit->first == -1) {
         while (lit != alnidx.end() && lit->first == -1) ++lit;
         seq1begin=lit->first; // should be 0
         seq2begin=lit->second;
         if (seq2begin == -1) {
            cerr << "Inside buildResult(), indel of the two aligned sequence follow each others, the algorithme might not be perfet for the linear space algorithm!\n";
            list<pair<int,int> >::const_iterator lix=lit;
            while (lix != alnidx.end() && lix->second == -1) {
               ++lix;
            }
            seq2begin=lix->second;
         }
      }
      else if (lit->second == -1) {
         while (lit != alnidx.end() && lit->second == -1) ++lit;
         seq1begin=lit->first;
         seq2begin=lit->second;  // should be zero
         if (seq1begin == -1) {
            cerr << "Inside buildResult(), indel of the two aligned sequence follow each others, the algorithme might not be perfet for the linear space algorithm!\nSeq1 start with -1\n";

            list<pair<int,int> >::const_iterator lix=lit;
            while (lix != alnidx.end() && lix->first == -1) {
               cerr << lix->first << "x" << lix->second << "  ";
               ++lix;
            }
            seq1begin=lix->first;
            cerr << "\nthe new seq1begin: " << seq1begin << endl;
         }
      }

      list<pair<int,int> >::reverse_iterator ri=alnidx.rbegin();
      if (ri->first == -1) {
         while (ri != alnidx.rend() && ri->first == -1) ++ri;
         seq1end=ri->first;
         seq2end=ri->second;
      }
      else if (ri->second == -1) {
         while (ri != alnidx.rend() && ri->second == -1) ++ri;
         seq1end=ri->first;
         seq2end=ri->second;
      }
   }
   buildAlnInfo(delta1, delta2);
   countnumgaps();
}



/* use two rows of array to memorize the coordinate
 * of the alignment starting point
 */
int LSDynaln::local() {
   alntype=LOCAL;
   allocmem();
   const int col=seq2->length()+1; // column of matrix
   const int row=seq1->length()+1;
   int i,j;
   // initialize match and indel score for row 0
   for (j=0; j< col; j++) {
      M[j]=IX[j]=IY[j]=0;
      // initialize position to default value -1
      MR[j<<1]=-1;  // faster version
      MR[(j<<1)+1]=-1;
   }
   M[col]=IX[col]=IY[col]=0;
   int *Mdia, *Mleft, *Mtop;
   int *IXdia, *IXleft, *IXtop;
   int *IYdia, *IYleft, *IYtop;
   // use MR array to remember the position
   // in the scoring matrix, the first column was set to -1

   int PLi, PLj, PCi, PCj; // position left and current
   //int *It;  // pointer to top insertion
   Smax=0;
   // for tracing the starting point of the alignment
   int *pos; // pointer to the position
   int Pmaxi, Pmaxj, score, scoreM, scoreIX, scoreIY;
   int top=0;
   int bottom=col;

   for (i=0; i < row-1; i++) {
      Mdia=M+top; Mleft=M+bottom;
      IXdia=IX+top; IXleft=IX+bottom;
      IYdia=IY+top; IYleft=IY+bottom;
      PLi=-1; PLj=-1;
      pos=MR+2;
      for (j=0; j < col-1; j++) {
         *(Mleft+1)=scoreM=max(*Mdia, max(*IXdia, *IYdia)) + ST.lookup(C1[i], C2[j]);
         *(IXleft+1)=scoreIX=max(*(Mdia+1)+gapi, *(IXdia+1)+gape);
         *(IYleft+1)=scoreIY=max(*Mleft+gapi, *IYleft+gape);
         score=max(max(scoreM, scoreIX), max(0,scoreIY));
         if (score > 0) {
            // not looking back on cell, the start position
            // is one before the actual alignment.
            if (score == scoreM) {
               PCi=*(pos-2); PCj=*(pos-1); // diagnal
            }
            else if (score == scoreIX) {
               PCi=*pos; PCj=*(pos+1);  // top
            }
            else if (score == scoreIY) {
               PCi=PLi; PCj=PLj; // left
            }
            else {
               cerr << "wrong condition in memorizing the starting point of local alignment\n";
               exit(1);
            }
            if (score > Smax) {
               Smax = score;
               Smaxi=i; Smaxj=j;
               Pmaxi=PCi; Pmaxj=PCj;
            }
         }
         else {
            PCi=i; PCj=j;
         }
         *(pos-2)=PLi; *(pos-1)=PLj; // left overwrite diag 
         PLi=PCi; PLj=PCj;       // current overwrite left
         pos += 2; // move two spaces to the right
         ++Mdia; ++Mleft;
         ++IXdia; ++IXleft;
         ++IYdia; ++IYleft;
      }
      M[top]=IX[top]=IY[top]=M[bottom]+gape;
      top=col-top;    // mechanism to reverse between 0 and 1
      bottom=col-bottom;
   }
   seq1begin=Pmaxi+1;
   seq1end=Smaxi;
   seq2begin=Pmaxj+1;
   seq2end=Smaxj;
   if (!alnidx.empty()) alnidx.clear();
   path(seq1begin, seq1end, seq2begin, seq2end);
   clearResult();
   return Smax;
}
LSDynaln::~LSDynaln() {
   //delete[] S; 
   //delete[] Itop;
   //delete[] M;
   //delete[] IX;
   //delete[] IY;
   // The base class members are done by base class
   delete[] MR;
   delete[] IXR;
   delete[] IYR;
}

// take top and bottom row pointers
void LSDynaln::debug_showmatrixForward(pair<int,int> &rows, const int col) {
   int *Matchptr, *IXptr, *IYptr;
   Matchptr = M + rows.first;
   IXptr = IX + rows.first;
   IYptr = IY + rows.first;
   int j;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
   Matchptr = M + rows.second;
   IXptr = IX + rows.second;
   IYptr = IY + rows.second;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
}
void LSDynaln::debug_showmatrixBackward(pair<int,int> &rows, const int col) {
   int *Matchptr, *IXptr, *IYptr;
   Matchptr = MR + rows.first;
   IXptr = IXR + rows.first;
   IYptr = IYR + rows.first;
   int j;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
   Matchptr = MR + rows.second;
   IXptr = IXR + rows.second;
   IYptr = IYR + rows.second;
   for (j=0; j<col; j++) {
      cout << Matchptr[j] << "," << IXptr[j] << "," << IYptr[j] << " | ";
   }
   cout << endl;
}
}
