#ifndef STDDEV_H
#define STDDEV_H

// (c) 1997 Kemin Zhou at The Molecular Sciences Institute

#include <utility>
#include <cmath>
#include <iostream>
#include <array>
#include <vector>

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
/**
 * Simple statistics for numbers
 */
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
      stddev(stddev&& other) = default;
      stddev& operator=(stddev&& other) = default;
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
       * read in n identical values of x
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
       * print avg,stddev,n to the output stream.
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
      double getTotal() const { return j*getMean(); }
		// postgres retuns the sample standard deviation
		pair<double, double> result() const { return make_pair(avg, sqrt(SS)); }
      /**
       * For tabular output.
       */
      friend ostream& operator<<(ostream &ous, const stddev &s) 
      {
         //ous << "mean\tstddev\tcount\n" << s.avg << '\t' << sqrt(s.var) << '\t' << s.j; 
         ous << s.avg << '\t' << sqrt(s.var) << '\t' << s.j; 
         return ous; 
      }

      /**
       * For human to read.
       */
      void show(ostream& ous) const {
         ous << avg << "+/-" << sqrt(var) << ":" << j;
      }

      /**
       * The object has any values accumulated.
       * @return true if this object has no data.
       */
      bool empty() const { return j==0; }
      /**
       * Rest this object to be used as a new object.
       */
      void clear() {
         j=0; SS=0; avg=0; var=0; 
      }
};

/**
 * Useful for working with multiple numbers.
 */
template<size_t N>
class Stddev {
	private:
		double var[N];
		double avg[N];
      /** the n value */
		int j;
		double SS[N];

	public:
      /**
       * Default constructor.
       */
		Stddev() : var{0}, avg{0}, j(0), SS{0} {}
      /**
       * Constructor from one element.
       */
		Stddev(double val[]) 
         : var{0}, avg{0}, j(1), SS{0} 
      {
         for (size_t i=0; i<N; ++i) avg[i] = val[i]; 
      }
      /**
       * Constructor from a std::array type.
       */
		Stddev(const array<double,N> val) 
         : var{0}, avg{0}, j(1), SS{0} 
      {
         for (size_t i=0; i<N; ++i) avg[i] = val[i]; 
      }
      /**
       * copy constructor. 
       * Since all types are primitive, no need for 
       * move constructor.  The default compiler generated
       * move constructor should work fine.
       */
      Stddev(const Stddev &other) 
         : var{0}, avg{0}, j(other.j), SS{0} 
      {
         for (size_t i=0; i<N; ++i) {
            var[i]=other.var[i];
            avg[i]=other.avg[i];
            SS[i]=other.var[i];
         }
      }

      Stddev(Stddev &&other) = default;
      /* wrong syntax
      Stddev(Stddev &&other) = default;
         : var(std::move(other.var)), 
            avg(std::move(other.avg)), 
            j(j), SS(std::move(other.SS))
      {
         other.j=0;
      }
      */

      Stddev& operator=(const Stddev& other) {
         if (this != &other) {
            for (size_t i=0; i<N; ++i) {
               var[i]=other.var[i];
               avg[i]=other.avg[i];
               j = other.j;
               SS[i]=other.var[i];
            }
         }
         return *this;
      }

      Stddev& operator=(Stddev&& other) = default;
      /*
      Stddev& operator=(Stddev&& other) {
         if (this != &other) {
            var = std::move(other.var);
            avg = std::move(other.avg);
            j = other.j;
            other.j = 0;
            SS = std::move(other.SS);
         }
         return *this;
      }
      */

      /**
       * To be used to accumulate input values.
       */
		void operator()(double x[]) {
			++j;
			if (j==1) { 
            for (size_t i=0; i<N; ++i) {
               avg[i] = x[i];
               var[i] = 0;
               SS[i] = 0;
            }
         }
			else { 
            for (size_t i=0; i<N; ++i) {
               double avgnew = avg[i] + (x[i]-avg[i])/j;
               SS[i] = (1-static_cast<double>(1)/(j-1))*SS[i] + j*pow(avgnew-avg[i], 2);
               var[i] = (j-1)*(var[i]/j + pow(avgnew-avg[i], 2));
               avg[i] = avgnew; 
            }
			}
		}

		void operator()(const array<double,N>& x) {
			++j;
			if (j==1) { 
            for (size_t i=0; i<N; ++i) {
               avg[i] = x[i];
               var[i] = 0;
               SS[i] = 0;
            }
         }
			else { 
            for (size_t i=0; i<N; ++i) {
               double avgnew = avg[i] + (x[i]-avg[i])/j;
               SS[i] = (1-static_cast<double>(1)/(j-1))*SS[i] + j*pow(avgnew-avg[i], 2);
               var[i] = (j-1)*(var[i]/j + pow(avgnew-avg[i], 2));
               avg[i] = avgnew; 
            }
			}
		}

      /**
       * Get mean value for the member with index i
       * use 0-based index
       */
      double getMean(size_t i) const { return avg[i]; }
      /**
       * @return the mean value.
       */
		const double* getMean() const { return avg; }
      array<double, N> getAverage() const {
         array<double, N> tmp;
         for (size_t i=0; i<N; ++i) {
            tmp[i] = avg[i];
         }
         return tmp;
      }
      double getStd(size_t i) const { return sqrt(var[i]); }
      /**
       * @return the standard deviation.
       */
		array<double, N> getStd() const { 
         array<double, N> tmp;
         for (size_t i=0; i<N; ++i) {
            tmp[i] = sqrt(var[i]);
         }
         return tmp;
      }

      double getSampleStd(size_t i) const { return sqrt(SS[i]); }
		array<double, N> getSampleStd() const { 
         array<double, N> tmp;
         for (size_t i=0; i<N; ++i) {
            tmp[i] = sqrt(SS[i]);
         }
         return tmp;
      }
      /**
       * Number of data points
       */
      int getCount() const { return j; }
		// postgres retuns the sample standard deviation
		pair<array<double,N>, array<double,N> > result() const { 
         return make_pair(getAverage(), getStd()); 
      }
      // better format, tabular format
      friend ostream& operator<<(ostream &ous, const Stddev<N> &s) 
      {
         //ous << "mean\tstddev\tcount\n" << s.avg << '\t' << sqrt(s.var) << '\t' << s.j; 
         for (size_t i=0; i<N; ++i) {
            ous << s.avg[i] << '\t' << sqrt(s.var[i]) << '\t';
         }
         ous << s.j; 
         return ous; 
      }
      /**
       * print count, [avg,stddev] ... for each feature
       * to the output stream.
       * will not output endl for flexibility.
       * @param delim delimter default is TAB.
       * @return ous output stream for chaining.
       */
      ostream& print(ostream &ous, const string &delim="\t") const {
         array<double, N> sd=getStd();
         ous << getCount();
         for (size_t i=0; i<N; ++i) {
            ous << delim << avg[i] << delim << sd[i];
         }
         return ous;
      }

      /**
       * The object has any values accumulated.
       */
      bool empty() const { return j==0; }
      /**
       * Reset this object to ground state.
       */
      void clear() {
         j=0; 
         for (size_t i=0; i<N; ++i) {
            SS[i]=0; avg[i]=0; var[i]=0; 
         }
      }
};

// use unsigned long
/**
 * The data must be two or more in size.
 * And must be sorted, the direction of the sorting is
 * unimportant.
 */
template<class T>
class Bisectsorted {
   private:
      vector<T> input;
      stddev avg1;
      stddev avg2;
      /**
       * pivot is the start of the uppper bound
       * [-----|P-----]
       *  123456
       *        789999
       */
      int pivot;

   public:
      Bisectsorted(const vector<T>& data) : input(data),
            avg1(input.front()), avg2(input.back()) { }
      Bisectsorted(vector<T>&& data) : input(std::move(data)),
            avg1(input.front()), avg2(input.back()) { }
      double getHighLowRatio() const {
         return avg2.getMean()/avg1.getMean();
      }
      double getLowHighRatio() const {
         return avg1.getMean()/avg2.getMean();
      }
      int getPivot() const {
         return pivot;
      }
      double getLog10Ratio() const {
         return log10(avg2.getMean()/avg1.getMean());
      }
      double getTotalHighLowRatio() const {
         return avg2.getTotal()/avg1.getTotal();
      }
      double getTotalLowHighRatio() const {
         return avg2.getTotal()/avg1.getTotal();
      }
      double getTotalLog10Ratio() const {
         return log10(getTotalHighLowRatio());
      }
      bool ratioGreater(double rcut) const {
         return getHighLowRatio() > rcut || getLowHighRatio() > rcut;
      }
      friend ostream& operator<<(ostream& ous, const Bisectsorted& bs) {
         ous << "[0:" << bs.input[0] << "--" << bs.pivot-1 << ":" << bs.input[bs.pivot-1] << "]["
            << bs.pivot << ":" << bs.input[bs.pivot] << "--" << bs.input.size()-1 << ":" 
            << bs.input.back() << "] low avg=" << bs.avg1 << " | high avg=" << bs.avg2 
            << " ratio: " << bs.getHighLowRatio()  << " " 
            << bs.getLowHighRatio() << endl;
         return ous;
      }
         
      /**
       * Separate numbers into two groups
       */
      void separate() {
         //  i      j
         //  ---=====
         int i=1;
         int j=input.size()-2;
         double diff1, diff2;
         while (j > i) {
            diff1 = abs(avg1.getMean() - input[i]);
            diff2 = abs(avg2.getMean() - input[i]);
            if (diff1 < diff2) {  // *--|
               avg1(input[i]);
               ++i;
               pivot=i;
            }
            else { // all numbers from here on belong to the larger group
               // this is the end of the algorithm
               pivot=i;
               while (i < j) {
                  avg2(input[i]);
                  ++i;
               }
            }
            if (i < j) {
               diff1 = abs(avg1.getMean() - input[j]);
               diff2 = abs(avg2.getMean() - input[j]);
               if (diff1 < diff2) {
                  // from start to j all belong to the lower group
                  pivot=j+1;
                  while (j > i) {
                     avg1(input[j]);
                     --j;
                  }
               }
               else {
                  pivot = j;
                  avg2(input[j]);
                  --j;
               }
            }
         }
         //cout << "Partition: [0:" << input[0] << "--" << pivot-1 << ":" << input[pivot-1] << "]["
         //   << pivot << ":" << input[pivot] << "--" << input.size()-1 << ":" << input.back()
         //   << " stat low: " << avg1 << " | stat high: " << avg2 << endl;
      }
};

}
#endif
