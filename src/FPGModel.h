#ifndef FPGMODEL_H
#define FPGMODEL_H

#include "GModel.h"

/** class: FPGModel
 * foot-print gene model 
 * is an extension of GModel and a specialization of the class
 * with fpkey.
 *
 * This class add a few parameter
 * from the footprint table that
 * records identity, coverage, and score
 *
 * This is used for Genewise pipeline and other
 * programs the consolidate blastx results into
 * concrete results.
 *
 */

namespace orpara {
class fpkey {
   public:
   fpkey(const string &q, const string &t, const int n) : qid(q), tid(t), fpnum(n) { }
   bool operator==(const fpkey &k) const { return qid==k.qid && tid==k.tid && fpnum==k.fpnum; }
   bool operator<(const fpkey &k) const {
      if (qid < k.qid) return true;
      if (qid > k.qid) return false;
      if (tid < k.tid) return true;
      if (tid > k.tid) return false;
      if (fpnum < k.fpnum) { return true; }
      return false;
   }
   bool operator>(const fpkey &k) const { return !(*this < k); }
   string qid, tid;
   int fpnum;
   friend ostream& operator<<(ostream &ous, const fpkey &k) {
      ous << k.qid << "\t" << k.tid << "\t" << k.fpnum;
      return ous; }
};

//template <class T> class GModel;
class FPGModel : public GModel<fpkey> {
   public:
      /* Take one row of footprint table, without quality
       * as input
       */
      FPGModel(const string &qid, const string &tid, int fpnum, 
            int ql, int tl, int qb, int qe, int tb, int te, 
            int nx, double ss, double ai, double qc, double  so)
         : GModel<fpkey>(fpkey(qid,tid,fpnum)), 
         qlen(ql), tlen(tl), qbegin(qb), qend(qe), tbegin(tb), tend(te),
         numexon(nx), sumscore(ss), avgiden(ai), qcov(qc), sumoverlap(so)
      { }
      friend ostream& operator<<(ostream& ous, const FPGModel &mm);
      /** 
       * check that all exons of mm is inside of exons of this object
       * Basic condition: mm inside this, mm score < this score
       * number of exon of this > number of exons of mm
       */
      bool cover(const FPGModel &mm) const;
      /**
       * ends similar, identity similar, score similar
       */
      bool similar(const FPGModel &mm, float fr=0.15) const;
      double getSumscore() const { return sumscore; }
      double getNormSumscore() const { return sumscore*(1-sumoverlap); }
      double getAvgiden() const { return avgiden; }

   private:
      int qlen, tlen, qbegin,qend,tbegin,tend,numexon;
      double sumscore,avgiden,qcov,sumoverlap;
};
}

#endif
