#include "DBStat.h"


Dbstat::Dbstat(const Dbstat &st) : numberdb(st.numberdb)  {
   stat = new SNavgstd[numberdb*numberdb];
   for (int i=0; i<numberdb*numberdb; i++) 
      stat[i]=st.stat[i];
}

Dbstat& Dbstat::operator=(const Dbstat &st) {
   if (this != &st) {
      if (numberdb<st.numberdb) {
         if (numberdb != 0) delete[] stat;
         numberdb=st.numberdb;
         stat=new SNavgstd[numberdb*numberdb];
      }
      for (int i=0; i<numberdb*numberdb; i++)
         stat[i]=st.stat[i];
   }
   return *this;
}
