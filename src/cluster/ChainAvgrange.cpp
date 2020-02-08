#include "ChainAvgrange.h"
#include "strformat.h"
#include <functional>
#include <iterator>

ChainAvgrange::~ChainAvgrange() {
	//cerr << "Calling ChainAvgrange destructor\n";
	for (unsigned int i=0; i<chain.size(); i++) {
		//cerr << "deleting chain element: " << *(chain[i]) << endl;
		delete chain[i];
	}
	//cerr << "all chain elements destroyed\n";
}

/** r will be saved into this object
 * This object will be responsible to delete the pointer.
 */
void ChainAvgrange::addRange(const alnrange* r) {
	if (chain.empty()) {
		chain.push_back( new avgrange(r) );
		return;
	}

	vector<avgrange*>::iterator vi, anchor, del;
	//vector<chain_iterator> needmerge;
	//for (int i=0; i<chain.size(); i++) {
	int nmcount=0;
	//for (vi=chain.begin(); vi != chain.end(); vi++) {
	vi=chain.begin();
	// merge while checking
	while (vi != chain.end()) {
		if ((*vi)->overlap(*r)) {
			++nmcount;
			if (nmcount==1) {
				anchor = vi;
				(*vi)->merge(r);
				vi++;
			}
			else if (nmcount > 1) {
				(*anchor)->merge(*vi);
				del=vi;
				vi++;
				delete *del;  // should we deallocate?
				chain.erase(del);
			}
			//mergedIterator.push_back(vi);
		}
		else vi++;
	}
	if (nmcount == 0) chain.push_back( new avgrange(r) );

	/* do the deletion after the checking id done
	if (mergedIterator.size() == 0) 
		chain.push_back( new avgrange(r) );
	else if (mergedIterator.size() > 1) {
		for (int i=1; i<mergedIterator.size(); i++) {
			(*mergedIterator[0])->merge(*mergedIterator[i]);
			chain.erase(mergedIterator[i]);
		}
	}
	*/
}

void ChainAvgrange::dbrows(list<string>& rows, const char dlm[]) const {
	for (size_t i=0; i<chain.size(); i++) {
		rows.push_back(seqid + dlm + itos(seqlen) + dlm
				+ chain[i]->asDelimitedString(dlm));
	}
}

ostream& operator<<(ostream &ous, const ChainAvgrange &ar) {
	for (unsigned int i=0; i<ar.chain.size(); i++) {
		ous << ar.seqid << "\t" << ar.seqlen << "\t";
		ar.chain[i]->writeTable(ous) << endl;
	}
	return ous;
}

int ChainAvgrange::checkSplits(list<string> &tablerows, ostream &ous) const {
	int numsplit=0;
	for (unsigned int i=0; i+1<chain.size(); i++) {
		SplitResult sres= chain[i]->testSplit(*chain[i+1]);
		if (sres.empty()) continue; 

		sres.setGuide(seqid,seqlen);
		// for database insertion
		list<string> rows=sres.outputRow(",");
		tablerows.insert(tablerows.end(), rows.begin(), rows.end());
		// for flat file dump
		rows=sres.outputRow("\t");
		copy(rows.begin(), rows.end(), ostream_iterator<string>(ous,"\n"));

		++numsplit;
	}
	return numsplit;
}

