#include "ncbitaxonomy.h"
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
   NCBITaxonomy taxtree;
   char *metagDataDir=getenv("METAG_DATA");
   string nodeFile = metagDataDir + string("/taxonomy/nodes.dmp");
   taxtree.buildTree(nodeFile);
   string nameFile = metagDataDir + string("/taxonomy/names.dmp");
   Taxon::loadTaxName(nameFile);
   string mappingFile = "germanyD_otu.tab";
   taxtree.countTaxon(taxtree.readMapping(mappingFile));
   vector<pair<string, int> > genusCount = taxtree.getGenusCount();
   cout << "genus hit_count\n";
   for (auto i=0; i<genusCount.size(); ++i) {
      cout << genusCount[i].first << " " << genusCount[i].second << endl;
   }
   vector<pair<string, int> > familyCount = taxtree.getFamilyCount();
   cout << "family hit_count\n";
   for (auto i=0; i<familyCount.size(); ++i) {
      cout << familyCount[i].first << " " << familyCount[i].second << endl;
   }

   return 0;
}

