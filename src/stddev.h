#ifndef STDDEV_H
#define STDDEV_H

// (c) 1997 Kemin Zhou at The Molecular Sciences Institute

#include <utility>
#include <cmath>
#include <iostream>

using namespace std;

/** 
 * use this module to get average and stddev from a set of values
 * This library is header only.
 * I used functional object.
 *
 * The formular that I used eliminated the under and over flow problem.
 * Average and Std can be calculted for any number.
 *
 * A straight implementation would have both under and over flow problems.
 */
namespace orpara {
class stddev {
	private:
		double var;  // variance
		double avg;
      /** the n value */
		int j;
		double SS;  // sample variance

	public:
      /**
       * Default constructor.
       */
		stddev() : var(0), avg(0), j(0), SS(0) {}
      /**
       * Constructor from one element.
       */
		stddev(double val) : var(0), avg(val), j(1), SS(0) {}
      /**
       * copy constructor. 
       * Since all types are primitive, no need for 
       * move constructor.  The default compiler generated
       * move constructor should work fine.
       */
      stddev(const stddev &other) : var(other.var), avg(other.avg), j(other.j),
         SS(other.SS) { }
      stddev& operator=(const stddev& other) {
         if (this != &other) {
            var=other.var;
            avg=other.avg;
            j=other.j;
            SS=other.SS;
         }
         return *this;
      }

      /**
       * To be used to accumulate input values.
       */
		void operator()(double x) {
			++j;
			if (j==1) { 
            avg = x; var = 0; SS=0; 
         }
			else { 
            double avgnew = avg + (x-avg)/j;
				SS = (1-static_cast<double>(1)/(j-1))*SS + j*pow(avgnew-avg, 2);
				var = (j-1)*(var/j + pow(avgnew-avg, 2));
				avg = avgnew; 
			}
			//cerr << "n: " << j << " avg =" << avg 
			//	<< " sqrt(SS)= " << sqrt(SS) << endl;
		}
      /**
       * read in n values of x
       */
      void operator()(double x, int n) {
         for (int i=0; i<n; ++i) {
            operator()(x);
         }
      }

      /**
       * @return the mean value.
       */
		double getMean() const { return avg; }
      /**
       * @return the standard deviation.
       */
		double getStd() const { return sqrt(var); }
		double getSampleStd() const { return sqrt(SS); }
      /**
       * print avg stddev n to the output stream.
       * will not output endl for flexibility.
       * @param delim delimter default is TAB.
       * @return output stream for chaining.
       */
      ostream& print(ostream &ous, const string &delim="\t") const {
         ous << avg << delim << getStd() << delim << j;
         return ous;
      }
      /**
       * Number of data points
       */
      int getCount() const { return j; }
		// postgres retuns the sample standard deviation
		pair<double, double> result() const { return make_pair(avg, sqrt(SS)); }
      friend ostream& operator<<(ostream &ous, const stddev &s) 
      {
         ous << "mean\tstddev\tcount\n" << s.avg << '\t' << sqrt(s.var) << '\t' << s.j; 
         return ous; 
      }
      /**
       * This object can be reused as empty after this operation.
       */
		void clear() {
         var=0; 
         avg=0;
         j=0;
         SS=0;
      }

      /**
       * The object has any values accumulated.
       */
      bool empty() const { return j==0; }
};
}
#endif
