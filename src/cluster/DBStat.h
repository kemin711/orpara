#ifndef DBStat_H
#define DBStat_H

#include <iostream>

using namespace std;

class SNavgstd {
   public:
   SNavgstd() { }
   SNavgstd(double sa, double ss, double na, double ns) 
      : scoreavg(sa), scorestd(ss), ngavg(na), ngstd(ns) { }
   double scoreavg, scorestd, ngavg, ngstd;
   friend ostream& operator<<(ostream &ous, const SNavgstd &sn) {
      ous << sn.scoreavg << " " << sn.scorestd << " " << sn.ngavg << " " << sn.ngstd; return ous; }

   double meanNgidentity() const { return ngavg; }
   double stddevNgidentity() const { return ngstd; }
   double stdNgidentity() const { return ngstd; }
};
/** This class stores the summary statistics of 
 * all x all (without self x self) SNavgstd elements.
 *
 * Internally I used a nxn array of SNavgstd to store the
 * information where n is the totall number of databases.
 */
class Dbstat {
   public:
      Dbstat() : numberdb(0) { }
      Dbstat(const int ndb) : numberdb(ndb) { 
         stat=new SNavgstd[numberdb*numberdb]; }
      Dbstat(const Dbstat &st);
      ~Dbstat() { delete[] stat; }
      /** index into db1*db2
       */
      void add(const int db1, const int db2, const SNavgstd &data)
      { stat[(db1-1)*numberdb + db2]=data; }

      Dbstat& operator=(const Dbstat &st);
      int getNumberOfDb() const { return numberdb; }
      /** return a pointer to the first element
       * db1_i => { db2_j | j=0, 1, ..., numdb-1 }
       *
       * You need to constrain yourself to getNumberOfDb()
       * elements of the array. The whole array is 
       * (numberdb * numberdb).
       */
      SNavgstd* getStatRow(const int db1) const {
         return stat + 1 + (db1-1)*numberdb; }

   private:

      int numberdb;
      /** it is numberdb x numberdb array of SNAvgstd
       * Not storing self x self comparisons.  In this
       * case you need to allocate (numberdb*numberdb)+1
       * elements.
       */
      SNavgstd* stat;
};

#endif
