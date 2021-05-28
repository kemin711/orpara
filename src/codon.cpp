#include "codon.h"
#include <cstring>
#include <cstdlib>
#include <fstream>

namespace orpara {

// this is the universal codon table
const string codon::univcodon = "ttt F ttc F tta L ttg L tct S tcc S tca S tcg S tat Y tac Y taa * tag * tgt C tgc C tga * tgg W ctt L ctc L cta L ctg L cct P ccc P cca P ccg P cat H cac H caa Q cag Q cgt R cgc R cga R cgg R att I atc I ata I atg M act T acc T aca T acg T aat N aac N aaa K aag K agt S agc S aga R agg R gtt V gtc V gta V gtg V gct A gcc A gca A gcg A gat D gac D gaa E gag E ggt G ggc G gga G ggg G";
//char codon::unknownaa='?';
char codon::unknownaa='X'; // this is the accepted letter 
// the following needs more generic treatment for portability
//char codon::codonfile[200]="/home/kzhou/etc/codontable.txt";
char codon::codonfile(string(DATADIR) + "/codontable.txt");

vector<map<string,char> > codon::codontables=vector<map<string,char> >(28);
vector<set<string> > codon::starts=vector<set<string> >(28);

// use a particular codong file
void codon::setCodonFile(const char file[]) {
   if (strlen(file) > 199) {
      cerr << "input path too long, should be less than 200 char\n";
   }
   //strcpy(codonfile, file);
   codonfile=file;
}

bool isCodonHeader(const string &line) {
   if (isdigit(line[0])) {
      if (line[1] == '.') return true;
      else {
         if (isdigit(line[1])) {
            if (line[2] == '.') return true;
            return false;
         }
      }
   }
   return false;
}
//void codon::readCodonTable(const string &file) {
void codon::readCodonTable() {
   cerr << "readCodonTable start\n";
   ifstream ins(codonfile);
   if (ins.fail()) {
      cerr << "Failed to open codon file: " << codonfile << endl;
      exit(1);
   }
   cout << "Reading codon tables from " << codonfile << endl;
   string line, startline;
   getline(ins, line);
   cout << line << endl << line[0] << " " << line[1] << endl;
   while (!ins.eof()) {
      if (line.empty()) getline(ins, line);
      else if (line[0]== '/' && line[1] == '*') {
         if (line[line.length()-2] == '*' && line[line.length()-1] == '/') {
            getline(ins, line);
         }
         else {
            getline(ins, line);
            while (true) {
               if (line.empty()) getline(ins,line);
               else if (line[0] == '*' && line[1] == '/') {
                  getline(ins,line);
                  break;
               }
               else getline(ins,line);
            }
         }
      }
      else if (line[0] == '/' && line[1] == '/') {
         getline(ins, line);
      }
      else if (isCodonHeader(line)) {
         string::size_type x=line.find('.');
         int codontable_id=atoi(line.substr(0,x).c_str());
         cout << "loading codon table " << codontable_id << endl;
         getline(ins,line);
         while (line.substr(0,5) != "  AAs") getline(ins,line);
         //line=line.substr(line.find('=')+2);
         line=line.substr(9);
         getline(ins, startline);
         startline=startline.substr(9);
         map<string,char> tmp;
         set<string> tmpstart;
         char base[4]={'T', 'C', 'A', 'G'};
         x=0;
         for (int i=0; i<4; ++i) {
            for (int j=0; j<4; ++j) {
               for (int k=0; k<4; ++k) {
                  string codonText=string(1,base[i]) + string(1, base[j]) + string(1,base[k]);
                  //cout << codonText << " => " << line[x] << endl;
                  //if (startline[x] == 'M' && codonText != "ATG") {
                  if (startline[x] == 'M') {
                     //cout << codonText << " is start\n";
                     tmpstart.insert(codonText);
                  }
                  tmp[codonText]=line[x++];
               }
            }
         }
         codontables[codontable_id]=tmp;
         starts[codontable_id]=tmpstart;
         getline(ins,line);
         while (!ins.eof() && !line.empty()) getline(ins,line);
      }
      else {
         cerr << line << endl << "shold not reach this condition for degug\n";
         exit(1);
      }
   }
      cerr << "readCodonTable done\n";
}
void codon::showAllCodonTables(ostream &ous) {
   if (codontables.empty()) readCodonTable();
   ous << "Valid codon tables: \n";
   for (unsigned int i=0; i<codontables.size(); ++i) {
      if (!codontables[i].empty()) 
         ous << i << "  ";
   }
   ous << endl;
}

//////////// end of global fuction ///////////////////


// given as a string of TTT F TTC F ....
// this function will parse this string
codon::codon(const std::string &def)
   : majorStart(), altstart(), tab() {
   strcpy(majorStart, "ATG");
   string::size_type i=0;
   while (i < def.size()) {
      tab[def.substr(i,3)] = def[i+4];
      i += 6;
   }
   convert();
}

codon::codon() 
   : majorStart(), altstart(), tab() {
   //static const string dft = "ttt F ttc F tta L ttg L tct S tcc S tca S tcg S tat Y tac Y taa X tag X tgt C tgc C tga X tgg W ctt L ctc L cta L ctg L cct P ccc P cca P ccg P cat H cac H caa Q cag Q cgt R cgc R cga R cgg R att I atc I ata I atg M act T acc T aca T acg T aat N aac N aaa K aag K agt S agc S aga R agg R gtt V gtc V gta V gtg V gct A gcc A gca A gcg A gat D gac D gaa E gag E ggt G ggc G gga G ggg G taa * tga * tag *";
   strcpy(majorStart, "ATG");
   unsigned int i=0;
   while (i<univcodon.size()) {
      tab[univcodon.substr(i,3)] = univcodon[i+4];
      i += 6;
   }
   convert();
   altstart.insert("ttg");
   altstart.insert("ctg");
   altstart.insert("atg");
}

// encode base into 2 bits 
int hashbase(const char n) {
   switch (n) {
      case 'a': case 'A':
      return 0;
      case 'c': case 'C':
      return 1;
      case 'g': case 'G':
      return 2;
      case 't': case 'T': case 'u': case 'U':
      return 3;
      case 'S': case 's':
      return 4;
      case 'W': case 'w':
      return 5;
      case 'R': case 'r':
      return 6;
      case 'Y': case 'y':
      return 7;
      case 'K': case 'k':
      return 8;
      case 'M': case 'm':
      return 9;
      case 'B': case 'b':
      return 10;
      case 'V': case 'v':
      return 11;
      case 'H': case 'h':
      return 12;
      case 'D': case 'd':
      return 13;
      case 'N': case 'n':
      return 14;
      default:
#ifdef DEBUG
      cerr << "Illegal base character: |" << n << "|(" << int(n) << ") inside hashbase\n";
#endif
      return 15;
   }
}

int hashcodon(const char c[3]) {
   if (strlen(c) == 1) return 64;
   else if (strlen(c) == 2) return 65;
   else {
      int c1=hashbase(c[0]);
      int c2=hashbase(c[1]);
      int c3=hashbase(c[2]);
      if (c1 > 3 || c2 > 3 || c3 > 3) {
         if ((c1 == 1 && c2 == 3) // Leu: CUN use CUA
             || (c1 == 2 && c2 == 3)  // Val: GUN
             || (c1 == 3 && c2 == 1)  // Ser UCN
             || (c1 == 1 && c2 == 1)  // Pro CCN
             || (c1 == 0 && c2 == 1)  // Thr ACN
             || (c1 == 2 && c2 == 1)  // Ala GCN
             || (c1 == 1 && c2 == 2)  // Arg CGN
             || (c1 == 2 && c2 == 2)  // Gly GGN
               ) { 
            c3=0;
         }
         else return 66;
      }
      return c1<<4 | c2 << 2 | c3;
   }
   //return hashbase( c[0] ) << 4  
    //     | hashbase( c[1] ) << 2
     //    | hashbase( c[2] ); 
}

int hashcodon(const string &cc) {
   if (cc.size() == 1) return 64;
   else if (cc.size() == 2) return 65;
   else {
      int c1=hashbase(cc[0]);
      int c2=hashbase(cc[1]);
      int c3=hashbase(cc[2]);
      if (c1 > 3 || c2 > 3 || c3 > 3) return 66;
      return c1<<4 | c2 << 2 | c3;
   }
}

// conver to the encoded codon table
void codon::convert() {
   map<string,char>::const_iterator it;
   for (it=tab.begin(); it != tab.end(); it++) {
      nuc2aa[hashcodon(it->first)]=it->second;
   }
   nuc2aa[64]='1';
   nuc2aa[65]='2';
   //nuc2aa[66]='?';
   nuc2aa[66]=unknownaa;
}

map<char, double> codon::getAAUniformFrequency() const {
   map<char,double> tmp;
   map<string,char>::const_iterator i;
   // first need to count the stop codons,
   // different codon tables may have different
   // number of stop codons!
   int numstops=0;
   for (i=tab.begin(); i != tab.end(); i++) {
      if (i->second == '*') ++numstops;
   }

   for (i=tab.begin(); i != tab.end(); i++) {
      if (i->second != '*') {
         tmp[i->second] += 1.0/(64.0-numstops);
      }
   }
   return tmp;
}

/* less efficient replace with more efficient implementation
char codon::operator[](const string &cd) {
   if (cd.size() > 3) {
      cerr << "codon should not be longer than 3\n";
      exit(1);
   }
   else if (cd.size() < 3) return ' ';

   string tmp;
   for (int i=0; i<cd.size(); i++) {
      tmp += tolower(cd[i]);
   }
   map<string,char>::const_iterator c=tab.find(tmp);
   if (c != tab.end()) return c->second;
   else return '?';
}
 */

/** more efficient version **/
char codon::operator[](const char cc[3]) {
   return nuc2aa[hashcodon(cc)];
}

char codon::operator[](const string &cd) {
   return nuc2aa[hashcodon(cd)];
}

void codon::use(const int tabid) {
   if (codontables[1].empty()) 
      readCodonTable();
   tab=codontables[tabid];
   if (tab.empty()) {
      cerr << "codon table " << tabid << " does not exist\n";
      showAllCodonTables(cerr);
      exit(1);
   }
   convert();
   altstart=starts[tabid];
}

void codon::show(ostream &ous) const {
   map<string,char>::const_iterator it;
   int i=1;
   for (it=tab.begin(); it != tab.end(); ++it) {
      ous << it->first << " " << it->second << "   ";
      if (i%4 == 0) ous << endl;
      ++i;
   }
   ous << "Major Start Codon: " << majorStart
      << " Alternative Start codons\n";
   set<string>::const_iterator s;
   for (s=altstart.begin(); s != altstart.end(); ++s) {
      ous << *s << " ";
   }
   ous << endl << endl;
}

}
