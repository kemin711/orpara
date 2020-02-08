#include <string>
#include <iostream>
#include <fstream>
#include <strformat.h>

using namespace std;

int main(int argc, char* argv[]) {
   string infile="/home/kzhou/work/chimera/MicpuN2/taxid_Mamiellaceae";
   ifstream IN(infile.c_str());
   if (IN.fail()) {
      cerr << "Failed to open " << infile << endl;
      exit(1);
   }
   string line, ids;
   getline(IN, line);
   vector<string> tmpids, allids;
   while (IN) {
      if (line[0] != '#' && line.size()>0) {
         //cout << line << endl;
         tmpids=dissect(line);
         if (tmpids.size()>0) {
            allids.insert(allids.end(), tmpids.begin(), tmpids.end());
         }
      }
      getline(IN, line);
   }
   vector<string>::const_iterator it=allids.begin();
   ids=*it;
   ++it;
   while (it != allids.end()) {
      ids += "," + *it;
      //cerr << *it << " ";
      ++it;
   }
   //cerr << endl;
//#IN >> ids;
 //  IN >> tmp;
  // while (!IN.eof()) {
   //   ids += ("," + tmp);
    //  IN >> tmp;
   //}
   //return ids;
   cerr << ids << endl;
   return 0;
}
