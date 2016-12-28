#ifndef MODEL_H
#define MODEL_H

#include <string>  
#include <vector> 
#include <list> 
#include <iostream>
                
using namespace std; 

namespace orpara {
class ExondirectionError : public exception {
   public:
      ExondirectionError() : message("Different direction of exon inside the same gene") { }
      ExondirectionError(const string &msg) : message(msg) { }
      const char* what() const throw() { return message.c_str(); }
      ~ExondirectionError() throw() { }
   private:
      string message;
};
             
/**          
 * A exon model of genes.
 * Only depends on standard library.
 */       
class Model { 
   public: 
      Model() : id(), exons() { }
      Model(const string &name) : id(name), exons() { }
      Model(const Model& other) : id(other.id), exons(other.exons) { }
      Model(Model &&other) : id(std::move(other.id)), exons(std::move(other.exons)) { }
      //void addExon(int bb, int ee) { exons.push_back(make_pair(bb,ee)); }
      void addExon(int bb, int ee) throw(ExondirectionError);
      bool operator==(const Model &mm) const;
      // return +, -, or ' ' if now known
      char direction() const;
      //pair<int,int> bound() const { return make_pair(exons[0].first, exons[exons.size()-1].second); }
      pair<int,int> bound() const { 
         return make_pair(exons.front().first, exons.back().second); }
      /* one model contains another model
       * Both ends ousdide of mm, and the intron of 
       * mm is a subset of this object.
       */
      bool contain(const Model &mm) const;
      /* only check for the introns, disregard the ends
       * This method is useful where some simple extension
       * of the partial models might have been wrong.
       */
      bool containIntron(const Model &mm) const;
      // checking only the begin and end of the model
      bool endsContain(const Model &mm) const;
      int numExons() const { return exons.size(); }
      //int start() const { return exons[0].first; }
      int start() const { return exons.front().first; }
      //int end() const { return exons[exons.size()-1].second; }
      int end() const { return exons.back().second; }
      // return all the intron boundaries
      vector<int> introns() const;
      friend ostream& operator<<(ostream& ous, const Model &mm);
      const string& getId() const { return id; }
 
   private: 
      string id; 
      //vector< pair<int,int> > exons;
      list< pair<int,int> > exons;
};  
}

#endif
