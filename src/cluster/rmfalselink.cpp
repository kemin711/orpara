#include <iostream>
#include <iomanip>
#include <fstream>
#include "mysql++.h"
#include <map>
#include "DBStat.h"
#include <gsl/gsl_cdf.h>
#include <cmath>
//#include <gsl/gsl_stats_double.h>


//#include "dbinfo.h"

/** given a set of organism genome databases (N)
 *  and its gene annotations.
 *  I did a all x all search in both directions.
 *  Then I picked the best one match from the 
 *  other databases. 
 *
 *  Query Database => (N-1) target databases.
 *
 *  There could be from 1 to N-1 possible hits.
 *  I am trying to use the global statistics to
 *  weed out week links.
 *
 *  In this case I have 16 genome. Only query
 *  with 5 ore more linkes are used to generate
 *  the statistics. Valance is defined as the
 *  number of targets found.  It is an indicator
 *  of how important this gene is relative to
 *  other organism. 
 *
 *  Vlance can be affected by random errors of 
 *  genome assembly and annotation. So we need
 *  a two dimontional system (valance, average
 *  identity) to classify a gene.
 *
 */

using namespace std;
using namespace mysqlpp;
/** This class is basically a convient and effient
 * implementation to weed out weak links
 * The private variables are used by the algorithm.
 */
class FitNgidentity {
   public:
      FitNgidentity() : numdb(0), mean(0), std(0), xsqltcut(0.09),
         avgmean(0), avgstd(0), scf(0), zarray(0), sczarray(0),
         scaledstdarray(0), stdslope(0.145), stdintercept(0.001)
         { }
      FitNgidentity(const int nd) 
         : numdb(nd), mean(0), std(0), xsqltcut(0.09),
            avgmean(0), avgstd(0), scf(0), zarray(0), sczarray(0),
            scaledstdarray(0), stdslope(0.145), stdintercept(0.001)
            { }
      //FitNgidentity(const int nd, const double *m, const double *s);
      ~FitNgidentity() { delete[] mean; delete[] std; delete[] scf; delete[] zarray; delete[] sczarray;  delete[] scaledstdarray; }
      /** given an array of SNavgstd it uses the ngidentity
       * part to initialize this object
       */
      void setRuler(const SNavgstd *snv);
      /** @param ngidens is an array of query to target
       * database ngidentity indexed by the (dbid-1)
       *  @param log is an output stream for log information.
       *  This is used for debuging phase only.
       * @return 1 for good, * 0 for gad, * -1 for modified lk: 
       *      removed some bad link -9 for cannot remove anymore.
       */
      int check(double *ngidens, ostream &log);
      double avgMean() const { return avgmean; }
      double avgStd() const { return avgstd; }
      void setCutoff(const double cf) { xsqltcut=cf; }
      /* floor(ngidentity*10)/2 hashing should be fine
      +-------+--------------------+
      | ngr20 | stddev(ngidentity) |
      +-------+--------------------+
      |     0 |         0.00549086 | should have not entry use 1
      |     1 |         0.04519686 |
      |     2 |         0.05633459 |
      |     3 |         0.05553033 |
      |     4 |         0.05011752 |
      |     5 |         0.00000000 | => 0.05011752 use 4
      +-------+--------------------+
      floor(ngidentity*10) 
      |     0 |         0.00549086 | [0.0, 0.1)
      |     1 |         0.02019686 | [0.1, 0.2)
      |     2 |         0.03519686 | [0.2, 0.3)
      |     3 |         0.04633459 | [0.3, 0.4)
      |     3 |         0.05553033 |
      |     4 |         0.05011752 |
      |     5 |         0.00000000 | => 0.05011752 use 4
      I should use y=0.1*x + 0.005
      */
      static double stdtable[6];


   private:
      int numdb;
      /** the index of the array is the dbid-1
       * We use 0-based index.
       * mean identity of db1 to all the db2's.
       */
      double *mean;
      /** array for the std for the whole population */
      double *std; 
      /** we should use this one for testing.  The whole population
       * is too broad. For every 20% range.
       * This is derived by scaling the stdtable.
       * scale factor = std[db2-1]/avgstd
       */
      double *scf; 
      /** for efficiency. allocated only once.*/
      double *zarray; 
      /** also for efficiency, stores scaled z-values. */
      double *sczarray; 
      /** for efficiency this should be computed only once.
       * */
      double *scaledstdarray;
      //double *adaptedstd;
      double xsqltcut;
      /** average over all means */
      double avgmean;
      double avgstd;
      double stdslope, stdintercept;
};

struct Resrow {
   Resrow() : db1(0), seq1(0), db2(0), ngiden(0) { }
   Resrow(int d1, unsigned int s1, int d2, double n) 
      : db1(d1), seq1(s1), db2(d2), ngiden(n) { }
   int db1;
   unsigned int seq1;
   int db2;
   double ngiden;
};

struct Intabcol {
   Intabcol() : db1id("db1id"), db2id("db2id"), 
      db1seqid("db1pid"), db2seqid("db2pid") { }
   string db1id, db2id, db1seqid, db2seqid;
};


void readStat(Connection &conn, const string &table, map<string, map<string, SNavgstd> > &dbindex);
//void readStatInt(Connection &conn, const string &table, SNavgstd* &stat, int &numdb);
/** use integer index instead of string map
 */
Dbstat* readStatInt(Connection &conn, const string &table);
//void check(const double* lk, const SNavgstd* qr, const int ndb);
int removeFalsePositiveHits(Connection &conn, const string &stattab, 
      Query &query, const string &outfile, vector<Resrow> &blk);
/**
 * The output table has the fixed schema:
db1 integer, seq1id integer unsigned, db2 integer, ngidentity float, 
primary key(db1,seq1id,db2))
*/
void storeBadInTable(Connection &conn, const vector<Resrow> &lk, const string &tab, const float lowestng);

int main(int argc, char* argv[]) {
   //MysqlDBInfo mydb=MysqlDBInfo::getAuthenInfo();
   //string stat_table="fungal5more_stat"; //round 1
   //string stat_table="fungalaln_stat2morehits"; // round 2
   string stat_table="fungal_db2dbstat_bidirection"; // round 3
   string database="kemin";
   string host="shake";
   string user="kzhou";
   string password="Hotcreek2800";
   string inputtable="fungalaln_top1_full";
   string badlinktab="badlinks";
   Intabcol icol;
   //Dbstat dbstat;
   int passed=0, failed=0, i;
   int trylimit=12;
   float lowngcut=0.44;

   i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-t")) trylimit=atoi(argv[++i]);
      else if (!strcmp(argv[i], "-s")) stat_table=argv[++i];
      else if (!strcmp(argv[i], "-i")) inputtable=argv[++i];
      else if (!strcmp(argv[i], "-o")) badlinktab=argv[++i];
      else if (!strcmp(argv[i], "-h")) host=argv[++i];
      else if (!strcmp(argv[i], "-d")) database=argv[++i];
      else if (!strcmp(argv[i], "--ming")) lowngcut=atof(argv[++i]);
      else {
      }
      ++i;
   }
   string outfile=badlinktab + ".tab";

   try {
      Connection conn(database.c_str(), host.c_str(), user.c_str(), password.c_str());
      Query query=conn.query();
      //query << "select db1, seq1name, db2, seq2name, score, ngidentity from "
       //  << inputtable << " order by db1, seq1name, db2 ";
       // Reverse direction clean up
      //query << "select db2, seq2name, db1, seq1name, max(score) as score, max(ngidentity) as ngidentity from "
       //  << inputtable 
        // << " group by db2, seq2name, db1 order by db2, seq2name, db1";
      query << "select " << icol.db1id << ","
         << icol.db1seqid << ", " << icol.db2id << ", "
         << icol.db2seqid << ", score, ngidentity"
         << " from " << inputtable 
         << " order by " << icol.db1id << ", "
         << icol.db1seqid << ", " << icol.db2id;
      vector<Resrow> badlinks;
      removeFalsePositiveHits(conn, stat_table, query, outfile, badlinks);
      storeBadInTable(conn,badlinks,badlinktab,lowngcut);
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      return 1;
   }

   return 0;
}
/** stat table should have the following definition:
 * mysql:kemin:shake>desc db2dbngstat95;
 * +-----------+------------------+------+-----+---------+-------+
 * | Field     | Type             | Null | Key | Default | Extra |
 * +-----------+------------------+------+-----+---------+-------+
 * | db1       | int(10) unsigned | NO   | PRI | NULL    |       |
 * | db2       | int(10) unsigned | NO   | PRI | NULL    |       |
 * | avgscore  | decimal(14,4)    | YES  |     | NULL    |       |
 * | stdscore  | double(26,4)     | YES  |     | NULL    |       |
 * | avgiden   | decimal(18,8)    | YES  |     | NULL    |       |
 * | stdiden   | double(34,8)     | YES  |     | NULL    |       |
 * | avgngiden | decimal(18,8)    | YES  |     | NULL    |       |
 * | stdngiden | double(34,8)     | YES  |     | NULL    |       |
 * | count     | bigint(21)       | NO   |     | 0       |       |
 * +-----------+------------------+------+-----+---------+-------+
 *
 * The sqlstr should pull the following columns
 * query << "select db1, seq1name, db2, seq2name, score, ngidentity from "
 *       << inputtable << " order by db1, seq1name, db2 ";
 *
 * This function will be able to work on both the forward and backward
 * direction.
 * The backward has multiple hits from the forward which we picked
 * only one hit. So the backward direction should use the max 
 * ngidentity.  We are elminating the hits if found lower that
 * other databases. Since we have done a 95% cutoff, this
 * should not be a problem.
 */
int removeFalsePositiveHits(Connection &conn, const string &stattab, 
      Query &query, const string &outfile, vector<Resrow> &blk) 
{
   Dbstat* dbstat;
   int passed=0, failed=0, i;
   int trylimit=12;

   ofstream OU(outfile.c_str());
   ofstream LOG("rmfalselink.log");

   try {
      dbstat=readStatInt(conn,stattab);
/** taken information from the raw data and do some test
*/
      double *links = new double[dbstat->getNumberOfDb()+1];
      double *linkscp = new double[dbstat->getNumberOfDb()+1];
      FitNgidentity fit(dbstat->getNumberOfDb());
      fit.setCutoff(0.09);
      StoreQueryResult res=query.store();
      StoreQueryResult::size_type r=0;
      res.disable_exceptions();
      //Row row=res.fetch_row();
      Row row=res[r];
      while (row) {
         int db1=row.at(0);
         /* db1_i (query) => {db2_1, db2_2, ..., db2_numdb} 
          * statistics (avg, std) */
         //SNavgstd* qrow = dbstat.getStatRow(db1);
         fit.setRuler(dbstat->getStatRow(db1));
         cout << "Query db " << db1 << endl;
         while (row && db1 == static_cast<int>(row.at(0))) {
            int seq1name=row.at(1);
            LOG << endl << db1 << " " << seq1name << ":\n----------\n";
            /* clear array */
            for (i=0; i<=dbstat->getNumberOfDb(); i++) {
               links[i]=0.0;
               linkscp[i]=0.0;
            }
            while (row && db1 == static_cast<int>(row.at(0)) 
                  && seq1name == static_cast<int>(row.at(1))) 
            {
               links[static_cast<int>(row.at(2))]=row.at(5); 
               linkscp[static_cast<int>(row.at(2))]=row.at(5); 
               //row=res.fetch_row();
               row=res[++r];
            }
            int rv = fit.check(links+1, LOG);
            int numtry=0;
            while (rv == -1) {
               if (numtry > trylimit) {
                  break;
               }
               ++numtry;
               LOG << " **** Round  #" << numtry << " Checking\n";
               rv=fit.check(links+1, LOG);
            }
            if (rv == 1) ++passed;
            else if (rv == 0) ++failed;
            else if (rv == -9) {
               LOG << "Cannot remove any more\n";
            }
            else { 
               LOG << " !!!! Number of try exceeded " << trylimit << " limit\n";
            }
            for (i=1; i<= dbstat->getNumberOfDb(); i++) {
               if (i == db1) continue;
               //if (linkscp[i] > 0.0) { 
               //   OU << db1 << "\t" << seq1name << "\t" << i << "\t" << linkscp[i] << endl;
               //}
               if (links[i] < 0.0) {
                  OU << db1 << "\t" << seq1name << "\t" << i 
                     << "\t" << linkscp[i] << "\n";
                  blk.push_back(Resrow(db1,seq1name,i,linkscp[i]));
               }
            }
         }
         LOG << endl;
      }
      delete[] links;
      delete[] linkscp;
      cerr << passed << " good matche groups q=>{t: 1-16}\n"
         << failed << " bad matches\n";
      delete dbstat;
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      return 1;
   }
   return 0;
}

void storeBadInTable(Connection &conn, const vector<Resrow> &lk, const string &tab, const float lowestng) {
   cerr << "storing " << lk.size() << " bad links ...\n";
   try {
      Query query=conn.query();
      query << "drop table if exists " << tab;
      query.execute();
      query.reset();
      query << "create table " << tab
         << "(db1 integer, seq1id integer unsigned, db2 integer, ngidentity float, primary key(db1,seq1id,db2))";
      query.execute();
      query.reset();
      query << "insert into " << tab << " values(%0, %1, %2, %3)";
      query.parse();
      for (int i=0; i<lk.size(); i++) {
         if (lk[i].ngiden >= lowestng) {
            query.execute(lk[i].db1, lk[i].seq1, lk[i].db2, lk[i].ngiden);
         }
      }
   }
   catch (exception &err) {
      cerr << err.what() << endl
         << " Failed to store bad links\n";
      exit(1);
   }
}


double FitNgidentity::stdtable[6]={ 0.03500686, 0.04219686, 0.05733459,
0.05653033, 0.05211752, 0.04811752 };


void FitNgidentity::setRuler(const SNavgstd *snv) {
   if (mean == 0) {
      mean=new double[numdb];
      std=new double[numdb];
      scf=new double[numdb];
      /* temporary for efficiency */
      zarray=new double[numdb];
      sczarray=new double[numdb];
      scaledstdarray=new double[numdb];
      //adaptedstd=new double[numdb];
   }
   double summ=0.0, sums=0.0;
   int i;
   for (i=0; i<numdb; i++) {
      mean[i]=snv[i].meanNgidentity();
      std[i]=snv[i].stdNgidentity();
      summ += mean[i];
      sums += std[i];
   }
   avgmean = summ / numdb; 
   avgstd = sums / numdb;
   for (i=0; i<numdb; i++) scf[i] = std[i]/avgstd;
   for (i=0; i<numdb; i++) {
      scaledstdarray[i]=stdslope*mean[i] + stdintercept;
      //cerr << "mean[i] " << mean[i] << " scaled std " << i << " " << scaledstdarray[i] << endl;
   }
}

int FitNgidentity::check(double *lk, ostream &log) {
   double chisqr=0.0;
   double zval, prob;
   //double *zarray= new double[numdb];
   double sumz=0.0, sumo=0.0, sumem=0.0, sumes=0.0;
   double minz=100.9, maxz=-100.1, maxng=-1.1, minng=2.2;
   //double *scaledstdarray = new double[numdb];
   int df=0;
   int i, maxngidx;
   log << "db_target\tzval\tobsngiden\tmeanngiden\tstdngiden\n";
   for (i=0; i<numdb; i++) {
      if (lk[i] > 0) {
         zval=(lk[i]-mean[i])/scaledstdarray[i];
         zarray[i]=zval;
         chisqr += zval*zval;
         ++df;
         sumz += zval;
         sumo += lk[i];
         sumem += mean[i];
         sumes += std[i];
         if (zval>maxz) maxz=zval;
         if (zval<minz) minz=zval;
         if (lk[i]>maxng) { maxng=lk[i]; maxngidx=i; }
         if (lk[i]<minng) { minng=lk[i]; }
         log << setw(3) << right << i+1 << "  "
            << setw(12) << left << showpos << setprecision(5) << zval 
            << setw(7) << left << noshowpos << setprecision(3) << lk[i] 
            << setw(6) << mean[i] << " +/- "
            << setprecision(4) << scaledstdarray[i] << endl;
      }
   }
   prob = gsl_cdf_chisq_P(chisqr, df);
   // discard single link if too bad with < 0.05 probility
   if (df == 1) {
      //if (sumz > -1.645) return 1;
      if (sumz > -1.5) return 1; // adaptive z
      if (maxng > 0.35) {
         double popz=(maxng-mean[maxngidx])/std[maxngidx];
         if (popz < -1.5) {
            log << "Removed, popzval=" << popz << endl;
            lk[maxngidx]=-1;
            return 0;
         }
         else return 1;
      }
      else {
         lk[maxngidx] = -1;
         return 0;
      }
   }

   log << "Chi Square=" << chisqr << " prob=" << prob << endl;
   if (prob <= xsqltcut) return 1;
   //(prob > xsqltcut) { scale then test
   log << " df=" << df << " average: zval=" << sumz/df
      << " ngiden=" << sumo/df
      << " \nexpected | " << sumem/df << " +/- " << sumes/df 
      << " | " << avgmean << " +/- " << avgstd 
      << endl;
   if (minz > -0.26 || (minz > -0.28 && minng > 0.36) || 
      (minz > -0.6 && minng > 0.375) ||
      (minz > -0.7 && minng > 0.395) ||
      (minz > -0.8 && minng > 0.52) ) {
      log << "Passed minz test\n";
      return 1;
   }
   // I am not sure we should discard low identities
   // they could be true for phylogenetic tree the 
   // low identities are reality for good measure
   if (maxz < -2 && maxng < 0.275) { 
      log << "all identity too low, discarded\n";
      for (i=0; i<numdb; i++) {
         if (lk[i] > 0.0000001) { lk[i] = -1; }
      }
      return 0;
   }
   // we need to scale both mean[i] and scaledstdarray[i]
   log << "Adapted Chi Square test\n";
   double sf = sumo/sumem; // scale ruler
   double stdadapt;
   chisqr=0.0;
   double scsumz=0.0, scminz=100.1, scmaxz=-100.1; 
   for (i=0; i<numdb; i++) {
      if (lk[i] <= 0.0000000001) continue;
      stdadapt=stdslope*mean[i]*sf + stdintercept;
      zval=(lk[i]-mean[i]*sf)/stdadapt;
      chisqr += zval*zval;
      sczarray[i]=zval;
      scsumz += zval;
      if (zval < scminz) scminz=zval;
      if (zval > scmaxz) scmaxz=zval;
      log << right << setw(3) << i+1 << "  "
         << setw(12) << left << showpos << setprecision(5) << zval 
         << setw(8) << left << noshowpos << setprecision(3) << lk[i] 
         << setw(8) << left << mean[i]*sf << " +/- "
         << setw (5) << setprecision(4) << left << stdadapt 
         << endl;
   }
   double scprob = gsl_cdf_chisq_P(chisqr, df);
   log << " Scaled prob=" << scprob << " avgz=" << scsumz/df << endl; 
   // I am going to remove by zval -1.65 will be removed
   int sctestrmcnt=0;
   for (i=0; i<numdb; i++) {
      if (lk[i] >= 0.0000000001 && (
               (zval < -1.65 && lk[i] < 0.375) ||
               (zval < -1.9 && lk[i] < 0.4) 
              )

      ) { // 0.05 probility
         lk[i]=-1;
         ++sctestrmcnt;
      }
   }
   if (sctestrmcnt > 0) return -1;
   if (scprob <= xsqltcut) { 
      log << "Passed scaled test with probability: " << scprob << endl;
      return 1; 
   }
   // try equal value
   double meanobsng=sumo/df;
   double tmpmeanstd=stdslope*meanobsng + stdintercept;
   // we want to have a tighter control
   tmpmeanstd *= 0.9;
   chisqr=0.0;
   log << "equal value test\n";
   for (i=0; i<numdb; i++) {
      if (lk[i] <= 0.0000001) continue;
      zval=(lk[i]-meanobsng)/(scf[i]*tmpmeanstd);
      chisqr += zval*zval;
      /*
      log << right << setw(3) << i+1 << "  " 
         << setprecision(5) << left << showpos << setw(8) << zval << "  " 
         << noshowpos << left << setprecision(3) << setw(6) << lk[i] << "  "
         << meanobsng << " +/- " << tmpmeanstd << endl;
         */
   }
   double evprob = gsl_cdf_chisq_P(chisqr, df);
   log << " Equal distance prob=" << evprob << endl;
   if (evprob <= xsqltcut) { 
      log << "PASSED equal distance test\n";
      return 1; 
   }

   // discard the worst low values and try again
   if (minz > -0.55)  {
      log << "all values seems to be high\n";
      return 1;
   }
   // the scaled maxz is more useful
   if (minng > 0.36 && (scmaxz > 1.24 || df < 3)) {
      log << "Passed: minimal ngidentity > 0.35 and scmaxz > 1.24\n";
      return 1;
   }
   if (scminz > -1 && minng > 0.35 && scprob < 0.22) 
      return 1;
   log << "trying to remove bad links ...\n";
   minz += 0.002;
   int removedCnt=0;
   /// remove according to first test
   for (i=0; i<numdb; i++) {
      if (lk[i] <= 0.0000000001) continue;
      if (maxz > 2.2 && zarray[i] < -1.5 && lk[i] < 0.75*maxng) {
         lk[i]=-1; ++removedCnt;
      }
      else if (zarray[i] < minz && lk[i] < 0.75*maxng && lk[i] < 0.49 && zarray[i] < -1.2) {
         lk[i]=-1; ++removedCnt;
      }
      else if ((zarray[i] < -0.5 || zarray[i] < minz) && lk[i] < 0.36 && lk[i] < 0.55*maxng) {
         lk[i]=-1; ++removedCnt;
      }
   }
   if (removedCnt > 0) return -1;

   /// use scaled test results
   scminz += 0.002;
   for (i=0; i<numdb; i++) {
      if (lk[i] <= 0.00000000001) continue;
      if (sczarray[i] < scminz) {
         if ((lk[i] < 0.75*maxng && lk[i] < 0.48)
               || lk[i] < 0.36) {
            lk[i]=-1; ++removedCnt;
         }
         else if (scminz < -1.65 && scmaxz < 1.6) {
            lk[i]=-1; ++removedCnt;
         }
      }
   }
   if (removedCnt == 0) { 
      log << "Cannot remove further\n";
      return -9;
   }
   return -1;
}

/** check the two arrays, ignore zero values in obsv 
 * obsv has n+1 elements, element[0] was not nused.
 * the index number is the dbid that start from 1.
 *
 * @param expt is an array containing two elements
 *    mean,std encoded in (2i-1, 2i), with i=1, 2, ..., 2n
 *    the array has 2n+1 elements, element at zero index
 *    was not used so that the index can start from 1.
 * */
double chisquare(const double *obsv, const double *expt, const int n, double* &zv) {
   double xsq = 0;
   double z;
   for (int i=1; i<=n; i++) {
      zv[i]=0; // clear old result
      if (obsv[i] > 0) {
         //z= (obsv[i] - exptmean[i])/exptstd[i];
         z= (obsv[i] - expt[2*i-1])/expt[2*i];
         zv[i]=z;
         xsq += z*z;
      }
   }
   return xsq;
}


/** the table structure for stat
+-----------+---------------+------+-----+---------+-------+
| Field     | Type          | Null | Key | Default | Extra |
+-----------+---------------+------+-----+---------+-------+
| db1       | varchar(30)   | NO   |     | NULL    |       |
| db2       | varchar(30)   | NO   |     | NULL    |       |
| avgscore  | decimal(14,4) | YES  |     | NULL    |       |
| stdscore  | double(26,4)  | YES  |     | NULL    |       |
| avgiden   | decimal(18,8) | YES  |     | NULL    |       |
| stdiden   | double(34,8)  | YES  |     | NULL    |       |
| avgngiden | decimal(18,8) | YES  |     | NULL    |       |
| stdngiden | double(34,8)  | YES  |     | NULL    |       |
+-----------+---------------+------+-----+---------+-------+

This is the string version. For better performance I am
implementing the array version.
*/

void readStat(Connection &conn, const string &table, map<string, map<string, SNavgstd> > &dbindex) {
   try {
      Query query=conn.query();
      query << "select db1, db2, avgscore, stdscore, avgngiden, stdngiden from "
         << table << " order by db1, db2";
      StoreQueryResult res=query.store();
      res.disable_exceptions();
      StoreQueryResult::size_type r=0;
      //Row row=res.fetch_row();
      Row row=res[r];
      //map<string, map<string,SNavgstd> > dbindex;
      int i=0;
      while (row) {
         string db1=string(row.at(0));
         map<string, SNavgstd> tmpmap;
         while (row && db1 == string(row.at(0))) {
            tmpmap.insert(make_pair(string(row.at(1)), SNavgstd(static_cast<double>(row.at(2)), static_cast<double>(row.at(3)), static_cast<double>(row.at(4)), static_cast<double>(row.at(5)))));
            //row=res.fetch_row();
            row=res[++r];
         }
         dbindex.insert(make_pair(db1, tmpmap));
      }
      cerr << dbindex.size() << " databases\n";
      map<string, map<string, SNavgstd> >::const_iterator it;
      map<string, SNavgstd>::const_iterator jt;
      for (it=dbindex.begin(); it != dbindex.end(); it++) {
         cerr << it->first  << "\n";
         for (jt=it->second.begin(); jt != it->second.end(); jt++) {
               cerr << jt->first << " " << jt->second << endl;
         }
      }
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
}

//void readStatInt(Connection &conn, const string &table, Dbstat &stat) {
Dbstat* readStatInt(Connection &conn, const string &table) {
   try {
      Query query=conn.query();
      query << "select count(distinct db1) from " << table;
      StoreQueryResult res=query.store();
      //Row row = res.fetch_row();
      StoreQueryResult::size_type r=0;
      Row row = res[r];
      //int totaldbs=row.at(0);
      int numdb =row.at(0);
      cerr << numdb << " databases\n";
      //SNavgstd* stat=new SNavgstd[totaldbs * totaldbs];
      //stat=new SNavgstd[numdb * numdb];
      Dbstat *stat = new Dbstat(numdb);

      query.reset();
      query << "select db1, db2, avgscore, stdscore, avgngiden, stdngiden from "
         << table << " order by db1, db2";
      res=query.store();
      res.disable_exceptions();
      r=0;
      //row=res.fetch_row();
      row=res[r];
      int i,j;
      while (row) {
         int db1=row.at(0);
         int db2=row.at(1);
         stat->add(db1,db2, SNavgstd(static_cast<double>(row.at(2)), static_cast<double>(row.at(3)), static_cast<double>(row.at(4)), static_cast<double>(row.at(5))));
         //row=res.fetch_row();
         row=res[++r];
      }
      return stat; // return local wrong
      /* debug show
      for (i=1; i<= totaldbs; i++) {
         for (j=1; j<= totaldbs; j++) {
            if (i != j) {
               cerr << i << " x " << j << endl;
               cerr << stat[i*j] << endl;
            }
         }
      }
      */
   }
   catch (exception &err) {
      cerr << err.what() << endl;
      exit(1);
   }
}
