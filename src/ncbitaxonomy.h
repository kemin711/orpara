#ifndef NCBITAXONOMY_H
#define NCBITAXONOMY_H

#include "gentree.h"
#include <string>
#include <map>
#include <vector>

using namespace std;

namespace orpara {
class Taxon {
   public:
      Taxon() : id(0), rank(-1), count(0) {}
      Taxon(int i) : id(i), rank(-1), count(0) {}
      /**
       * Constructor
       * @param i taxid
       * @param r rankid
       */
      Taxon(int i, int r) : id(i), rank(r), count(0) {}
      Taxon(int i, int r, int c) : id(i), rank(r), count(c) {}
      /**
       * Constructor
       * @param i taxid
       * @param rn rank name in string format. Typical values are
       *     species, genus, family, superkindom, ...
       */
      Taxon(int i, string rn) : id(i), rank(getRankId(rn)), count(0) {}
      int getId() const { return id; }
      string& getName() const { return id2scientific[id]; }
      /**
       * increment the count field.
       */
      void addCount(int cnt) { count += cnt; }
      int getCount() const { return count; }
      int getRank() const { return rank; }
      bool rankHigherThan(const Taxon &other) { 
         return taxRankOrder[getRankName(rank)] < taxRankOrder[getRankName(other.getRank())]; }
      /**
       * to see the rank represented by the rank number in this
       * object is higher in the order than the rank number given.
       */
      bool rankHigherThan(int rnum) { 
         return taxRankOrder[getRankName(rank)] < taxRankOrder[getRankName(rnum)]; }
      bool rankHigherThan(const string &rname) { 
         return taxRankOrder[getRankName(rank)] < taxRankOrder[rname]; }
      bool rankLowerThan(const Taxon &other) { 
         return taxRankOrder[getRankName(rank)] > taxRankOrder[getRankName(other.getRank())]; }
      bool rankLowerThan(int rnum) { 
         return taxRankOrder[getRankName(rank)] > taxRankOrder[getRankName(rnum)]; }
      bool rankLowerThan(const string &rname) { 
         return taxRankOrder[getRankName(rank)] > taxRankOrder[rname]; }
      bool sameRank(const Taxon &other) { 
         return rank == other.rank; }
      bool sameRank(const string &rname) { 
         return rank == getRankId(rname); }
      bool sameRank(int rnum) { 
         return rank == rnum; }

      /**
       * Tax ran in sorted order, for binary look up.
       * The Taxon will use the integer version for saving 
       * memory.
       * "no rank" will be represented by -1. This is a special
       * value.
       */
      static int getRankId(const string &rn);
      static string getRankName(int rankid) { if (rankid==-1) { return "no rank"; } else { return taxRanks[rankid]; } }
      /** 
       * Loading only scientific name from the names.dmp file
       */
      static void loadTaxName(const string &infile);
      /**
       * taxid starts from 1. This is the node.dmp file
       * @return the taxid of the taxname according to the 
       *     NCBI taxonomy database. notaxid (-1) for not found. 
       * May have to switch to unsigned long int
       */
      static int getTaxid(const string &taxname);
      static const int notaxid = -1;

   private:
      /**
       * Taxid of this taxon in the NCBI taxonomy
       */
      int id;
      int rank;
      int count;
      
      /**
       * Helper function for getRankId()
       */
      static int binary_search_rank(const string& item, int left, int right);
      /**
       * "no rank" will not be part of the vector.
       */
      static vector<string> taxRanks;
      /**
       * Look up table to order the rank from high (samller number)
       * to low (larger number). This can be used for tree
       * traversal up to certain levels.
       */
      static map<string, int> taxRankOrder;
      /**
       * Scientific name of taxon to id mapping.
       * For look up.
       */
      static map<string, int> scientific2id;
      /**
       * Id to scientific name look up
       */
      static map<int, string> id2scientific;
};

class TaxonVisitor {
   public:
      TaxonVisitor() : result() { }
      void operator()(node<Taxon>* n) {
         if (n->data.getCount()>0) result.push_back(n->data);
      }
      vector<Taxon> getResult() const { return result; }
      vector<pair<string,int> > getTaxonCount() const;

   private:
      vector<Taxon> result;
};


/**
 * This should be mostly a singleton class
 * TODO: write a function to searilize the taxonomy tree
 * so that the loading can be very fast instead of paring
 * the text input from NCBI.
 */
class NCBITaxonomy : GenericTree<Taxon> {
   public:
      NCBITaxonomy() : GenericTree(), indexById() { }
      /**
       * Build the internal tree from input file
       * @param infile is the NCBI dump file
       */
      void buildTree(const string &infile);
      //void countTaxon(vector<pair<string, int> > &taxCount);
      /**
       * Count taxons if given a vector or list
       * of tuples<string, int, string>
       * @param taxcnt abundance of each taxon given by the last element.
       *   The input is (cluster_id, depth, species+strain). Sometimes
       *   the species+strain may be higher or lower than species.
       *   The taxon description may be ambiguous.
       *
       * This is the input function.
       */
      void countTaxon(const vector<tuple<string, int, string> > &taxcnt);
      void updateCount(int taxid, int cnt);
      vector<pair<string, int> > getRankCount(const string &rk) const;
      vector<pair<string, int> > getGenusCount() const { return getRankCount("genus"); }
      vector<pair<string, int> > getFamilyCount() const { return getRankCount("family"); }

      static string removeSquare(const string &str);
      /**
       * command line or debug function for testing the class
       * Read OTU taxonomy mapping from the 16S pipeline and
       * @param file is the otu mapping file with 15 columns
       * @return the result as a vector.
       * The plotter could also implement this.
       *
       * No side-effect on this object.
       * The result is used by countTaxon as input
       */
      static vector<tuple<string, int, string> > readMapping(const string &file);

   private:
      /**
       * For quickly locate the Taxon node by tax id in the 
       * taxonomy tree.
       */
      map<int, node<Taxon>* > indexById;
};
}
#endif
