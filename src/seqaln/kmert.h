#ifndef KMER_H
#define KMER_H

#include <string>
#include <vector>
#include <iostream>
#include <exception>
#include <cmath>
#include <bitset> // for debug print out
#include <iterator>
#include <cassert>
//#include <fstream>
#include <stdexcept>
#include "kmerhelper.h"
#include <set>
#include <map>
#include <algorithm>
// why kerm.h?
#include "kmer.h"
#include "kmerbase.h"

using namespace std;

namespace orpara {
/**
 * Template version of kmer.  Kmer length
 * is a parameter in this class. K should be 
 * less than 16.
 *
 * This class is for individual sequences.
 * This will be dedicated to DNA 4 bases
 * K is a Nontype templae parameter.
 * This class emphasize the location information.
 * Kmercount class emphasizes the count without
 * storing the location information.
 * @see KmerCount
 */
template<unsigned int K> 
class Kmert : public KmerBase<K> {
   private:
      /**
       * Each instance has one mask based on K parameter. 
       * Each K class share one mask
       */
      //static unsigned int mask;

      /**
       * Not sure we should store the original sequence
       */
      string seq;
      /**
       * hashed kmer interger value for position from 0 to L-k+1
       * This is tied to the input sequence length.
       */
      vector<uint64_t> hashval;
      /** 
       * The locations of each kmer in the input sequence.
       * First index is the kmer's hash value.
       * Vector store the 0 or more location of the kmer's 
       * first base's 0-based index in seq. Empty vector
       * signifies the lack of kmer in the entire sequence.
       * Each index from 0 to loc.size() is the hash value of
       * each possible kmer.
       */
      vector<vector<unsigned int> > loc;
      /**
       * location of the reverse complement
       */
      vector<vector<unsigned int> > locrc;

      /**
       * DNA fragements after cuting the DNA with
       * each kmer. There are multiple kmers.
       * Each kemer is indexd by their hash value.
       */
      mutable vector<vector<int> > frag;
      mutable vector<vector<int> > fragrc;
      mutable bool digested;

      /**
       * Compute the kmer for all positions
       * from 0 to seq.length()-k+1
       */
      void kmer2int();
      void tabulateLoc();
      void init();

   public:

      /**
       * hashval is initialized with empty vector.
       */
      Kmert() : seq(), hashval(), loc(mask+1), locrc(mask+1), 
         frag(loc.size()), fragrc(loc.size()), digested(false) {}
      /**
       * hashval is intialized with -1.
       * Pow(4,k) is the same as 2<<(2*k-1)
       * There are pow(4,k) possible k-mers with ACGT at symbol.
       * Each K-mer is represented with a integer number
       * from 0 to pow(4,k)-1
       */
      Kmert(const string &s)
         : seq(s), hashval(s.length()-K+1, -1), 
            loc(mask+1), locrc(mask+1), 
            frag(loc.size()), fragrc(loc.size()), digested(false)
      {
         init();
      }

      /**
       * fill the frag vector. Break up the input
       * sequence into fragments.
       * This operation should be done only on request.
       */
      void digest() const;

      /** for testing */
      void showLocation() const;
      /** testing function
       */
      void showFragment() const;
      /**
       * Debug function
       */
      void showHashval() const {
         copy(hashval.begin(), hashval.end(), ostream_iterator<int>(cout, ", "));
         cout << endl;
      }

      /**
       * obtain the result
       */
      const vector<uint64_t>& getHashValue() const {
         return hashval;
      }

      bool isPalindrome(pair<int,int> &lpreg, vector<double> &freq, vector<double> &rcfreq) const;

      /*
      unsigned int revcompBits(unsigned int hv) {
         //cout << bitset<32>(hv) << " " << hv << " input\n";
         unsigned int hvc = (~hv);
         //cout << bitset<32>(hvc) << " " << hvc << " complement\n";
         size_t i;
         unsigned int res=(3&hvc);
         for (size_t i=0; i<K-1; ++i) {
            res <<= 2;
            hvc >>= 2;
            res |= (3&hvc);
         }
         //cout << bitset<32>(res) << " " << res << " after reverse\n";
         //cout << bitset<32>(res&mask) << " " << (res&mask) << " after clean up\n";
         return res&mask;
      }
      */

      int getSequenceLength() const { return seq.length(); }
      const string& getSequence() const { return seq; }
      /**
       * Christopher added this one.
       * Not sure seq is from this object or external.
       * If external this function should be static.
      vector<double> cal_KmerFrequency(const string &seq);
      This one does not even compile, so leave out of the 
      compilation stage.
       */

      /**
       * convert base to integer
       */
      /**
       * Set the mer size (K) and build the mask
       * It create a bit string with 2*K 1 from the right.
       *   |-> 2*K  <-|
       * |0|1|1|....|1|
       * set mask to 2*K 1's
       */
      //static void buildMask();
};


//template<int K>
//unsigned int Kmert<K>::mask = computeMask(K);

///////// template kmer class implementation ////////////

/*
template<int K>
void Kmert<K>::buildMask() {
   mask=(2<<(2*K-1)) - 1;
   // for debug
   cout << "mask value: " << bitset<32>(mask) << endl;
}
*/

template<int K>
void Kmert<K>::kmer2int() {
   // calculate kmer integer hash value for all kmers in the sequence
   size_t i;
   unsigned int v=0;
   // build the intial hash value for the first kmer.
   for (i=0; i<K-1; ++i) {
      v <<= 2;
      v |= base2int(seq[i]);
   }
   for (i=0; i<seq.length()-K+1; ++i) {
      v<<=2;
      v |= base2int(seq[i+K-1]);
      v &= mask;
      hashval[i]=v;
   }
}

template<int K>
void Kmert<K>::tabulateLoc() {
   for (size_t i=0; i<hashval.size(); ++i) {
      if (hashval[i] != -1) {
         loc[hashval[i]].push_back(i);
      }
   }
   //vector<int> tmp(hashval.size());
   //cout << "inside tabulateLoc revcomp hash values\n";
   //for (unsigned int i=hashval.size()-1; i != 0; --i) {
   //   cout << revcompBits(hashval[i]) << ", ";
  // }
   //cout << revcompBits(hashval[0]) << endl;

   // convert to rc position
   //cerr << "loc size: " << loc.size() << " locrc size: "
   //   << locrc.size() << endl;
   //   https://www.biostars.org/p/113640/ check this one
   for (unsigned int i=mask; int(i) > -1; --i) {
      if (!loc[i].empty()) {
         for (unsigned int j=loc[i].size()-1; int(j)> -1; --j) {
            //locrc[revcompBits(i)].push_back(seq.length()-K+1-loc[i][j]);
            locrc[revcompKmerInt(i)].push_back(seq.length()-K+1-loc[i][j]);
         }
      }
   }
}

template<int K> void Kmert<K>::init() {
   kmer2int();
   tabulateLoc();
   //digest();
}

template<int K>
void Kmert<K>::showLocation() const {
   // i is the kmer hash value
   cout << "location on the top strand\n";
   for (size_t i=0; i<loc.size(); ++i) {
      cout << i << ":  ";
      if (!loc[i].empty()) 
         copy(loc[i].begin(), loc[i].end(), ostream_iterator<int>(cout, ", "));
      cout << endl;
   }
   cout << "location on the reverse strand\n";
   for (size_t i=0; i<locrc.size(); ++i) {
      cout << i << ":  ";
      if (!locrc[i].empty()) 
         copy(locrc[i].begin(), locrc[i].end(), ostream_iterator<int>(cout, ", "));
      cout << endl;
   }
   cout << endl;
}

template<int K>
void Kmert<K>::digest() const {
   if (digested) return;
   for (size_t i=0; i<loc.size(); ++i) {
      int previous=0;
      frag[i].clear();
      for (size_t j=0; j<loc[i].size(); ++j) {
         frag[i].push_back(int(loc[i][j])-previous);
         previous=loc[i][j];
      }
      frag[i].push_back(int(seq.length())-previous);
   }
   for (size_t i=0; i<locrc.size(); ++i) {
      int previous=0;
      fragrc[i].clear();
      for (size_t j=0; j<locrc[i].size(); ++j) {
         fragrc[i].push_back(int(locrc[i][j])-previous);
         previous=locrc[i][j];
      }
      fragrc[i].push_back(int(seq.length())-previous);
   }
   digested = true;
}

template<int K>
void Kmert<K>::showFragment() const {
   if (!digested) digest();
   for (size_t i=0; i<frag.size(); ++i) {
      cout << i << '\t';
      copy(frag[i].begin(), frag[i].end(), ostream_iterator<int>(cout, "\t"));
      cout << endl;
   }
   cout << "fragments on the reverse complement\n";
   for (size_t i=0; i<fragrc.size(); ++i) {
      cout << i << '\t';
      copy(fragrc[i].begin(), fragrc[i].end(), ostream_iterator<int>(cout, "\t"));
      cout << endl;
   }
}

/* not even compile, wrong usage of KmerCount class
template<int K>
vector<double> Kmert<K>::cal_KmerFrequency(const string &seq) {
   KmerCount mycounter(5);
   mycounter(seq, 1);
   vector<double> kmer_freq=mycounter.getFrequency();
   return kmer_freq;
}

template<int K>
double Kmert<K>::calculateDistance(vector<double>& freq1, vector<double>& freq2){
   
    double sum=0;
    for(int i=0;i<freq1.size();++i){

    	 if(freq1[i]>0){

     		 double diff = pow(freq1[i]-freq2[i], 2);
      		 sum += diff;  
       				 
     			

         }

 
    }
    return sum;
}
*/

template<int K>
bool Kmert<K>::isPalindrome(pair<int,int> &lpreg, vector<double> &freq,vector<double> &rcfreq) const {
   // kmer -> number of occurance
   vector<int> kmerCount(loc.size());
   int both=0, single1=0, single2=0;
   //cerr << loc.size() << " kmers\n";
   //cout << "\nsingle on top strand\n";
   set<int> loop;
   set<int>  stem;
   vector<unsigned int> myloop;
   std::vector<unsigned int>::iterator it;
   int cc=0;
   int count=0; 
   int count1=0;
   for (size_t i = 0; i<loc.size(); ++i) {
        if (loc[i].size() > 0 && locrc[i].size() > 0) {
         kmerCount[i] = 2;
         ++both;
  	        
         if(loc[i].size()==locrc[i].size())
             stem.insert(loc[i].begin(), loc[i].end());

      
         }
        else if (loc[i].empty() && !locrc[i].empty()) {
           ++single2;
           kmerCount[i] = 1;
        }
        else if (!loc[i].empty() && locrc[i].empty()) {
            ++single1;
            kmerCount[i] = 1;
            loop.insert(loc[i].begin(), loc[i].end());
      }
      else {
         kmerCount[i] = 0;
      }
   }
   float dr = float(both)/(both + single1);

  // cout<<dr<<endl;
      if (dr > 0.6) {
           
        vector<pair<int,int> > stemreg = combineRanges(stem,2);
       	
        vector<pair<int,int> > stemcombinereg;
     	vector<pair<int,int> > stemcombinereg_ten;

     	vector<pair<int,int> > left_stem;
     	vector<pair<int,int> > right_stem;


     	pair<int, int> newpair;
     	vector<pair<int,int>>::iterator it = stemreg.begin();
     	vector<int> stems;

     	std::map <int,int> stem_pos;

     	int count_ten=0;
     
     if(stemreg.size()<=1)
       return false;
      while( it != stemreg.end())
      {
   
    		 
  	         newpair.first  =(*it).first;
  	         newpair.second =(*it).second;
                 ++it;
                 stems.push_back(newpair.second-newpair.first);

                 if(newpair.second-newpair.first>=10)
  			stemcombinereg_ten.push_back(newpair);
                 stem_pos[stems.size()-1]=newpair.first;
		

      }
    if(stemcombinereg_ten.size()<=1)
         return false;
   
     
     int left=0;
     int right=stemcombinereg_ten.size()-1;

     int totall=0;
     
 
     while(left<right){

           int l1=stemcombinereg_ten[left].second-stemcombinereg_ten[left].first+1;

           int l2=stemcombinereg_ten[right].second-stemcombinereg_ten[right].first+1;

  
           if(abs(l1-l2)<4){
      //      cout<<stemcombinereg_ten[left].first<<","<<stemcombinereg_ten[left].second<<" "<<stemcombinereg_ten[right].first<<","<<stemcombinereg_ten[right].second<<endl;

           if(stemcombinereg_ten[left].second<stemcombinereg_ten[right].first){
            
                pair<int,int> left_pair,right_pair;

                left_pair.first= stemcombinereg_ten[left].first;
	        left_pair.second= stemcombinereg_ten[left].second;
                right_pair.first= stemcombinereg_ten[right].first;
                right_pair.second= stemcombinereg_ten[right].second;
                left_stem.push_back(left_pair);
                right_stem.push_back(right_pair);
 

           }
          
               ++left;
	       --right;
               ++totall;

           }


          else if(l1>l2+2){
  
              int c=0;
         
              while(c<3){
                 
                  --right;
                  int l3=stemcombinereg_ten[right].second-stemcombinereg_ten[right].first+1;
                 // cout<<stemcombinereg_ten[right].first<<stemcombinereg_ten[right].second<<endl;
                  if(abs(l1-l3)>1){
 			++c;
                        if(c==3){
			  ++left;
                          right=right+c-1;
                        }
                        
                  }
                  else{
                    // cout<<stemcombinereg_ten[left].first<<","<<stemcombinereg_ten[left].second<<" "<<stemcombinereg_ten[right].first<<","<<stemcombinereg_ten[right].second<<endl;
                     if(stemcombinereg_ten[left].second<stemcombinereg_ten[right].first){
            
             		pair<int,int> left_pair,right_pair;

            	        left_pair.first= stemcombinereg_ten[left].first;
	    	        left_pair.second= stemcombinereg_ten[left].second;
                        right_pair.first= stemcombinereg_ten[right].first;
                        right_pair.second= stemcombinereg_ten[right].second;
             	        left_stem.push_back(left_pair);
             	  	right_stem.push_back(right_pair);
 

          	    }


                       ++left;
		       --right;
                       ++totall;
 			break;

	         }
             
              }

          }


          else if(l2>l1+2){

		int c=0;
                while(c<3){

                 ++left;
                 int l4=stemcombinereg_ten[left].second-stemcombinereg_ten[left].first+1;
                 if(abs(l2-l4)>1){
 			++c;

                        if(c==3){
			  --right;
			  left=left-c+1;
                        }
                        
                  }
                  else{
                //   cout<<stemcombinereg_ten[left].first<<","<<stemcombinereg_ten[left].second<<" "<<stemcombinereg_ten[right].first<<","<<stemcombinereg_ten[right].second<<endl;

		 if(stemcombinereg_ten[left].second<stemcombinereg_ten[right].first){
            
             		pair<int,int> left_pair,right_pair;
             		left_pair.first= stemcombinereg_ten[left].first;
	    	        left_pair.second= stemcombinereg_ten[left].second;
             		right_pair.first= stemcombinereg_ten[right].first;
             		right_pair.second= stemcombinereg_ten[right].second;
            	        left_stem.push_back(left_pair);
             		right_stem.push_back(right_pair);
 
 
          	 }


                       ++left;
		       --right;
                       ++totall;
 			break;
                    

	         }               

          }

    } 
}
 
   int lbegin = left_stem[left_stem.size()-1].first;

   int lend = right_stem[right_stem.size()-1].second;

   std::set<int>::iterator it1;

   set<int> myloop;
   
   for (it1 =loop.begin(); it1 != loop.end(); it1++) {
     
     if(*it1>=lbegin && *it1<=lend){

        myloop.insert(*it1);

        loop.erase(it1);


     }

  
   }

   
    vector<pair<int,int> > loopreg = combineRanges(loop,2);
    vector<pair<int,int> > loopreg2 = combineRanges(myloop,2);

    

   if(loopreg2.size()>0){

       cout<<"The loop starts at "<<loopreg2[0].first<<","<<" and ends at "<<loopreg2[loopreg2.size()-1].second<<endl;


       int loop_size= loopreg2[loopreg2.size()-1].second-loopreg2[0].first+1;

       /* Not even compile, added by Christopher commenting it out
       if(loop_size>50){

             int l1=loopreg2[0].first-1;
             int l2= loopreg2[loopreg2.size()-1].second+1;

             vector<double> freql=cal_KmerFrequency(seq.substr(0,l1+1);
             vector<double> freqm=cal_KmerFrequency(seq.substr(l1+1,l2-l1-1));
             vector<double> freqr=cal_KmerFrequency(seq.substr(l2,seq.length()-1-l2));

              double sum1= calculateDistance(freql,freq);
              double sum2= calculateDistance(freql,rcfreq);

              double sum3= calculateDistance(freqm,freq);
              double sum4= calculateDistance(freqm,rcfreq);

             
   	      double sum5= calculateDistance(freqr,freq);
              double sum6= calculateDistance(freqr,rcfreq);

              
	      if(((sum1<sum2) && (sum5<sum6)) || ((sum1>sum2) && (sum5>sum6))){
                   
                   return false;
              
              }
                   
              else{

                   


	

              }             


       }
       */
   }
   else {
      cout<<"The loop size is 0 and the sequence is a perfect palindrome"<<endl;
      return true;
   }   
      return true;
   }
   return false;
}
}
#endif
