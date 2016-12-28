#include <iostream>
#include <string>
#include <vector>
#include <codon.h>
#include <bioseq.h>

using namespace std;

int main(int argc, char* argv[]) {
   cout << "testing linker problem\n";

   string str("ksksksk");
   cout << str << endl;
   str = "string from literial";
   cout << str << endl;
   vector<string> vs;
   vs.push_back("ksksksskskksskks Frist");
   vs.push_back("Second String after Frist");
   vs.push_back("Third string no room");
   vector<string>::size_type i;
   for (i=0; i<vs.size(); ++i) {
      cout << vs[i] << "  ";
   }
   cout << endl;
   cout << "done\n";

   bioseq seq("sometitle", "TAGTGATCCCCCCCCCCCCCJJJJJJJLLLLLLLLKKKKACGTGAWWWWEWAAAACCCVVVNNNMMMKKKLLL");
   cout << seq << endl;
   codon cod;
   cout << cod["AAA"] << cod["ATG"] << endl;

   return 0;
}
