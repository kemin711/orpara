#ifndef CHAINAVGRANGE_H
#define CHAINAVGRANGE_H

#include "alnrange.h"
#include <vector>
#include <iostream>

/*
 * To buiild a chain of average ranges for chimera detection
 */

using namespace std;

class ChainAvgrange {
	public:
		ChainAvgrange(const string &id, int len) 
			: seqid(id), seqlen(len), chain() {}
      /* this object will deallocate the memory given by the
       * addRange method
       */
		~ChainAvgrange();

		/* either r merge with one of the chain element
		 * or create a new avgrange object out of r
		 */
		void addRange(const alnrange* r);
		int size() const { return chain.size(); }

		/* fill the list with rows of db formated
		 * quoted string values.
		 */
		void dbrows(list<string>& rows, const char dlm[]=",") const;
		/* tab-dlimited output, good for database loading */
		friend ostream& operator<<(ostream &ous, const ChainAvgrange &ar);
		bool isChimera() const { return chain.size()>1; }

		/* Accumulate results to the tablerows, also output the
		 * result in tab-delimited format to the output stream
		 */
		int checkSplits(list<string> &tablerows, ostream &ous) const;
		typedef vector<avgrange*>::iterator chain_iterator;

	private:
		string seqid;
		int seqlen;
		vector<avgrange*> chain;
};

#endif
