#include "ncbitaxonomy.h"
#include <fstream>
#include <iostream>
#include <strformat.h>
#include <list>
#include <iterator>

namespace orpara {
// in sorted order for binary look up
vector<string> Taxon::taxRanks={
   "class", "family", "forma", "genus", "infraclass", "infraorder",
   "kingdom", "order", "parvorder", "phylum", "root", "species", 
   "species group", "species subgroup", "subclass", "subfamily", 
   "subgenus", "subkingdom", "suborder", "subphylum", "subspecies", 
   "subtribe", "superclass", "superfamily", "superkingdom", 
   "superorder", "superphylum", "tribe", "varietas"
};

map<string, int> Taxon::taxRankOrder = {
   {"root", 0},
   {"superkingdom", 1}, {"kingdom", 2}, {"subkingdom", 3},
   {"superphylum", 4}, {"phylum", 5}, {"subphylum", 6},
   {"superclass", 7}, {"class", 8}, {"subclass", 9}, {"infraclass", 10}, 
   {"superorder", 11}, {"order", 12}, {"parvorder", 13}, {"suborder", 14}, {"infraorder", 15}, 
   {"superfamily", 16}, {"family", 17}, {"subfamily", 18}, 
   {"tribe", 19}, {"subtribe", 20}, 
   {"genus", 21}, {"subgenus", 22}, 
   {"species group", 23}, {"species subgroup", 24}, 
   {"species", 25}, {"subspecies", 26}, 
   {"varietas", 27}, {"forma", 28}, {"norank", 99}};

int Taxon::binary_search_rank(const string &item, int left, int right) {
   if (left > right) return -1;
   int middle = (left + right)/2;
   if (taxRanks[middle] == item) 
      return middle;
   if (item < taxRanks[middle]) {
      return binary_search_rank(item, left, middle-1);
   }
   return binary_search_rank(item, middle+1, right);
}

int Taxon::getRankId(const string &rn) {
   if (rn == "no rank") return -1;
   int idx = binary_search_rank(rn, 0, taxRanks.size()-1);
   if (idx == -1) {
      throw runtime_error("rank: " + rn + " is unknown");
   }
   return idx;
}

map<string, int> Taxon::scientific2id={};
map<int, string> Taxon::id2scientific={};

void Taxon::loadTaxName(const string &infile) {
   ifstream inf(infile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + infile);
   }
   map<string,int>::iterator n2id_it;
   map<int,string>::iterator id2n_it;
   map<string,int> uniqueName2id;
   map<int,string> id2uniqueName;
   vector<tuple<int, string, string> > withUnique;
   set<string> duplicatedName; // duplicated name
   set<int> duplicatedId;

   // need to do a 2 pass to get the duplicated entry 
   // replaced
   string line, name;
   getline(inf, line);
   int cnt=0;
   while (!inf.eof()) {
      if (line.substr(line.length()-2) == "\t|") {
         line.resize(line.size()-2);
      }
      vector<string> row = split(line, "\t|\t");
      ++cnt;
      if (row[3] == "scientific name") {
         name = row[1];
         if (!row[2].empty()) {
            withUnique.push_back(make_tuple(stoi(row[0]), row[1], row[2]));
         }
         id2n_it = id2scientific.find(stoi(row[0]));
         if (id2n_it == id2scientific.end()) {
            id2scientific[stoi(row[0])]=name;
         }
         else { // ids are unique, you should not get here
            cout << "id: " << row[0] << " has multiple scientific names: "
               << id2n_it->second << " and " << name << "\n";
         }
         n2id_it = scientific2id.find(name);
         if (n2id_it == scientific2id.end()) {
            scientific2id[name]=stoi(row[0]);
         }
         else {
            if (row[2].empty()) {
               //throw runtime_error("multiple name: " + name + " taxid: " + row[0] + " had no unique name ");
               //cerr << "Warning: multiple name: " << name << " taxid: " << row[0] << " had no unique name\n";
               // make up unique name by add an unique integer
               withUnique.push_back(make_tuple(stoi(row[0]), row[1], row[1] + " " + row[0]));
            }
            duplicatedName.insert(name);
            duplicatedId.insert(n2id_it->second);
            //cout << "scientific name: " << name << " has multiple ids: "
            //   << n2id_it->second << " and " << row[0] << "\n";
         }
      }
      getline(inf, line);
   }
   cerr << cnt << " rows of data\n";
   // second pass to deal with duplicated names
   // 1 remove duplicated name 
   for (auto x=duplicatedName.begin(); x != duplicatedName.end(); ++x) {
      scientific2id.erase(*x);
   }
   // 2. remove duplicated id
   for (auto x=duplicatedId.begin(); x != duplicatedId.end(); ++x) {
      id2scientific.erase(*x);
   }
   // 3. add duplicated using unqiue names rather than names
   cout << withUnique.size() << " has unique names\n";
   cnt=0;
   for (auto i=0; i<withUnique.size(); ++i) {
      if (duplicatedName.find(get<1>(withUnique[i])) != duplicatedName.end()) {
         //cout << "adding unique name " << get<0>(withUnique[i]) << " "
         //   << get<1>(withUnique[i]) << " " << get<2>(withUnique[i]) << endl;
         ++cnt;
         scientific2id.insert(make_pair(get<2>(withUnique[i]), get<0>(withUnique[i])));
         id2scientific.insert(make_pair(get<0>(withUnique[i]), get<2>(withUnique[i])));
      }
   }
   cerr << cnt << " unique names added " << scientific2id.size() << " scientific names "
      << id2scientific.size() << " tax ids loaded from " << infile << "\n";
}

int Taxon::getTaxid(const string &taxname) {
   map<string,int>::const_iterator it = scientific2id.find(taxname);
   if (it == scientific2id.end()) {
      return notaxid;
   }
   else {
      return it->second;
   }
}

////////////// NCBITaxonomy  Class ////////////////////

string NCBITaxonomy::removeSquare(const string &str) {
   //cout << str << " => " << deleteChr(deleteChr(str, '['), ']') << endl;
   return deleteChr(deleteChr(str, '['), ']');
}

void NCBITaxonomy::buildTree(const string &infile) {
   cerr << "Building tree from " << infile << endl;
   ifstream inf(infile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + infile);
   }
   string line;
   //map<int, node<Taxon>* > indexById;
   //map<int, pair<int, string> > pending; // cannot construct now
   // id, parent, rank_id
   list<tuple<int, int, int> > pending; // cannot construct now
   getline(inf, line); // first line is speciel and is the root
   vector<string> row = split(line, "\t|\t");
   if (row[0] != "1") {
      throw runtime_error("first line must be root");
   }
   root->setData(Taxon(1));
   indexById[1]=root;
   // first round will only make the node whose parent is 1
   getline(inf, line);
   map<int, node<Taxon>* >::iterator it;
   node<Taxon>* ptr;
   int numNodes=1;
   while (!inf.eof()) {
      if (line.substr(line.length()-2) == "\t|") {
         line.resize(line.length()-2);
      }
      vector<string> row = split(line, "\t|\t");
      it = indexById.find(stoi(row[1]));
      if (it != indexById.end()) {
         //cerr << row[0] << " " << row[1] << " " << row[2] << " is indexById " << endl;
         ptr = it->second;
         if (ptr->child == 0) {
            ptr->addChild(Taxon(stoi(row[0]), Taxon::getRankId(row[2])));
            indexById.insert(make_pair(ptr->child->data.getId(), ptr->child));
         }
         else {
            node<Taxon>* newNode = ptr->child->addSibling(Taxon(stoi(row[0]), Taxon::getRankId(row[2])));
            indexById.insert(make_pair(newNode->data.getId(), newNode));
         }
      }
      else {
         pending.push_back(make_tuple(stoi(row[0]), stoi(row[1]), Taxon::getRankId(row[2])));
      }
      getline(inf, line);
      ++numNodes;
   }
   cerr << numNodes << " total nodes. First pass done " << indexById.size() << " nodes indexById\n";
   // now finish the rest
   list<tuple<int,int,int> >::iterator lit, x;
   int round=1;
   while (!pending.empty()) { // finish all the rest
      cout << "round " << round << ": " << pending.size() << " pending, "
         << indexById.size() << " indexById\n";
      lit=pending.begin();
      while (lit != pending.end()) {
         //cout << "working on node: " << get<0>(*lit) << " " << get<1>(*lit) << " " 
         //   << get<2>(*lit) << endl;
         it = indexById.find(get<1>(*lit));
         if (it != indexById.end()) {
            ptr = it->second; // parent node found
            if (ptr->child == 0) {
               ptr->addChild(Taxon(Taxon(get<0>(*lit), get<2>(*lit))));
               //cout << "add to indexById\n";
               indexById.insert(make_pair(get<0>(*lit), ptr->child));
            }
            else {
               node<Taxon>* newNode = ptr->child->addSibling(Taxon(get<0>(*lit), get<2>(*lit)));
               indexById.insert(make_pair(get<0>(*lit), newNode));
            }
            x=lit;
            ++lit;
            pending.erase(x);
         }
         else {
            //cout << "node: " << get<0>(*lit) << " is waiting for its parent node: " 
            //   << get<1>(*lit) << " to be indexById\n";
            ++lit;
         }
      }
      ++round;
   }
   cerr << "tree construction done\n";
}

void NCBITaxonomy::updateCount(int taxid, int cnt) {
   node<Taxon>* p = indexById[taxid];
   p->data.addCount(cnt);
   p = p->parent;
   while (p != root) {
      p->data.addCount(cnt);
      p = p->parent;
   }
   root->data.addCount(cnt);
   //cerr << root->data.getCount() << " total count in the tree\n";
}

void NCBITaxonomy::countTaxon(const vector<tuple<string, int, string> > &taxcnt) {
   cerr << "counting taxon hits from " << taxcnt.size() << " nodes\n";
   string taxname;
   for (auto i=0; i<taxcnt.size(); ++i) {
      taxname = get<2>(taxcnt[i]);
      //cerr << "working on " << i << " " << taxname << endl;
      //string taxname = removeSquare(taxname); // some taxname may have square bracket, 
      // we are removing it here. This could be done by the supplier (caller).
      int taxid;
      if ((taxid=Taxon::getTaxid(taxname)) != Taxon::notaxid
            || (taxid=Taxon::getTaxid(taxname=chopend(taxname, " uncultured"))) != Taxon::notaxid
            || (taxid=Taxon::getTaxid(taxname=chopend(taxname, "uncultured"))) != Taxon::notaxid 
            || (taxid=Taxon::getTaxid(taxname=chopLastWord(taxname))) != Taxon::notaxid // one extra chopping
            || (taxid=Taxon::getTaxid(taxname=chopLastWord(taxname))) != Taxon::notaxid
            || (taxid=Taxon::getTaxid(taxname=firstword(taxname))) != Taxon::notaxid) 
      {
         //cout << "updating taxid " << taxid << " " << taxname << endl;
         updateCount(taxid, get<1>(taxcnt[i]));
      }
      else {
         cerr << "Warning: " << get<0>(taxcnt[i]) << " " << taxname << " has no taxid\n";
      }
   }
}

vector<tuple<string, int, string> > NCBITaxonomy::readMapping(const string &file) {
   cerr << "Reading OTU taxon mapping from " << file << endl;
   ifstream ifs(file);
   if (ifs.fail()) throw runtime_error("cannot fine file: " + file);
   string line, clname; 
   getline(ifs, line); // header
   getline(ifs, line); // first data line
   // clid, depth, species+strain
   vector<tuple<string, int, string> > om;
   vector<string> row, lastrow;
   while (!ifs.eof()) {
      //cout << line << endl;
      row = split(line, "\t");
      string species=row[4];
      if (row[5] != "\\N") {
         species += " " + row[4];
      }
      species = removeSquare(species);
      if (!lastrow.empty() && row[0] == lastrow[0]) {
         if (stoi(row[2]) >= stoi(lastrow[2])) { // replace last one with this one
            om.back() = make_tuple(row[0], stoi(row[2]), species);
         }
         else { // saving the last one
            cout << "multiple mapping for cluster " << row[0] << endl
               << "ignoring the weaker mapping\n" << line << endl
               << "saving\n";
            copy(lastrow.begin(), lastrow.end(), ostream_iterator<string>(cout, " | "));
            cout << endl;
         }
      }
      else {
         om.push_back(make_tuple(row[0], stoi(row[2]), species));
      }
      getline(ifs, line);
      lastrow = row;
   }
   cout << om.size() << " OTUs\n";
   // there is no need to sort the taxons
   return om;
}

///////////// visitor class ////////////////////////////////////////
vector<pair<string,int> > TaxonVisitor::getTaxonCount() const {
   vector<pair<string,int> > tmp;
   for (auto i=0; i<result.size(); ++i) {
      tmp.push_back(make_pair(result[i].getName(), result[i].getCount()));
   }
   return tmp;
}

////////////////////////////////////////////////////////////////////////

void visitRank(node<Taxon> *n, TaxonVisitor &vis, const string &r) {
   if (n != 0) {
      if (!n->isRoot() && n->data.rankLowerThan(r)) return; // not vising lower ranks
      if (n->data.sameRank(r)) {
         vis(n);
      }
      visitRank(n->child, vis, r);
      visitRank(n->sibling, vis, r);
   }
}

vector<pair<string, int> > NCBITaxonomy::getRankCount(const string &rk) const {
   TaxonVisitor visitor;
   visitRank(root, visitor, rk);
   return visitor.getTaxonCount();
}
}
