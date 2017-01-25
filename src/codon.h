#ifndef CODON_H
#define CODON_H

// (c) 2002 Kemin Zhou at orapra.com

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <set>

using namespace std;

namespace orpara {
/** convert base to integer
 * A=>0, C=>1, G=>2, T,U=>3, 
 * S,W,R,Y,K,M,B,V,H,D,N => 4-14
 * regardless of case
 * For efficient coding decoding
 */
int hashbase(const char n);
/** 
 * convert 3-letter codon to a integer
 * from 0-63 for regular codon with letter from [ACGT]
 * if codon length is 1 then return 64, 2 ten 65, other 66
 */
int hashcodon(const char c[3]);
int hashcodon(const string &cc);

/** 
 * Used for translation of nucleic acids. Normally we use T for U. 
 * We are kind of doing DNA CDS, not directly translating mRNA.
 */
class codon {
   public:
      /** initialize with the default universal codon table. */
      codon();
      /** use a given codon table 
       * given as a string of TTT F TTC F ....
       * This class should be able to input from more 
       * formats, such as files etc.
       * @param def is the string version of the codontable. 
       *   It is in 3nt codon followed by space followed by sigle letter amino
       *   acid code.
       * */
      codon(const std::string &def);
      /** switch to a particular codon table.
       * for example useTable(12). This will switch to Candida yeast
       * codon table. */
      void use(const int tabid);

      /** translate codon into one letter aa code
       * This is the key method to convert codon into aa.
       */
      char operator[](const string &cd);

      /** translate codon into amino acid one-ltter code.
       * This is the efficient version.
       * unknown codons such as NNN will be translated into '?'
       * This is controlled by the unknowaa variable.
       * 
       * Use hashcodon helper function. If codon is 1 nt,
       * the codon integer 64 is '1', 2 nt return 2
       * for codons containing unknow bases, it return
       * 'unknown aa' that is controlled by the static variable:
       *    unknownaa.
       */
      char operator[](const char cc[3]);

      map<char,double> getAAUniformFrequency() const;
      /** universal codon
       * It was initialized with lower case 3 letter codons followed by single
       * letter amino acids code.  Total 64 codons with star represents stop
       * codon.
       **/
      static const string univcodon;
      /** Symbol to use for unknow amino acid default X.
       * translation of codons with Ns, such as NNN NNA, .... Even two
       * bases determine an AA, we still translate it to X.  
       **/
      static char unknownaa;
      //static void readCodonTable(const string &file);
      static void readCodonTable();
      static void setCodonFile(const char file[]);
      /** this is a global parameter it should be assigned to the system-wide
       * directory. Now default to user or programmer's home directory
       * $HOME/etc/codontable.txt
       */
      static char codonfile[200]; // $HOME/etc/codontable.txt
      /** for debug purpose */
      void show(ostream &ous) const;
      static void showAllCodonTables(ostream &ous);

   private:
      char majorStart[4]; //"ATG";
      /** alternative starts */
      set<string> altstart;
      /** 
       * temporary used by constructors. 
       * String foramt => tab ==convert()==>array format [nuc2aa]
       * This is the current codon table used.
       * String version.
       * 3-letter CODON => One letter AA Code
       */
      map<string, char> tab;
      /**
       * Is a one-dimensional array of single upper case letters.
       * use hashcodon to index the codons. This is for fast look up.
       * The single letter is determined by the codon table 
       * map<string, char> tab used.
       */
      char nuc2aa[67]; // [64] 1, [65] 2 [66] last character is ?
      /** 
       * convert the encoded table into one dimentional array 
       * for fast look up.
       * unknown codon use ?, 1 for 1 Nucleotide at end of seq.
       * 2 for two nucleotide of codon left even 6 aa can
       * be determined by 2-nucleotides.
       */
      void convert();

      /** 
       * so far there are 23 with 7,8 removed
       * actually 21 table. We will leave 7,8 empty.
       * I gave 28 for future growth.
       * To see this table is loaded, you need to test
       * the codon table at index 1.
       */
      static vector<map<string,char> > codontables;
      static vector<set<string> > starts;
};

// the universal codon table is defined like this
// char nuc2aa[64]={ /* AAA */  'K', /* AAC */  'N', /* AAG */  'K',
//   /* AAU */  'N', /* ACA */  'T', /* ACC */  'T', /* ACG */  'T',
//   /* ACU */  'T', /* AGA */  'R', /* AGC */  'S', /* AGG */  'R',
//   /* AGU */  'S', /* AUA */  'I', /* AUC */  'I', /* AUG */  'M',
//   /* AUU */  'I', /* CAA */  'Q', /* CAC */  'H', /* CAG */  'Q',
//   /* CAU */  'H', /* CCA */  'P', /* CCC */  'P', /* CCG */  'P',
//   /* CCU */  'P', /* CGA */  'R', /* CGC */  'R', /* CGG */  'R',
//   /* CGU */  'R', /* CUA */  'L', /* CUC */  'L', /* CUG */  'L',
//   /* CUU */  'L', /* GAA */  'E', /* GAC */  'D', /* GAG */  'E',
//   /* GAU */  'D', /* GCA */  'A', /* GCC */  'A', /* GCG */  'A',
//   /* GCU */  'A', /* GGA */  'G', /* GGC */  'G', /* GGG */  'G',
//   /* GGU */  'G', /* GUA */  'V', /* GUC */  'V', /* GUG */  'V',
//   /* GUU */  'V', /* UAA */  '*', /* UAC */  'Y', /* UAG */  '*',
//   /* UAU */  'Y', /* UCA */  'S', /* UCC */  'S', /* UCG */  'S',
//   /* UCU */  'S', /* UGA */  '*', /* UGC */  'C', /* UGG */  'W',
//   /* UGU */  'C', /* UUA */  'L', /* UUC */  'F', /* UUG */  'L', 
//   /* UUU */ 'F'};
// Any exception must be handled outside this table.
}
#endif
