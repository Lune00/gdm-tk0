#ifndef _network_hpp
#define _network_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include "sample.hpp"
#include "inter2d.hpp"
#include "groupRelationData.hpp"
#include "fsafe.hpp"
#include "near.hpp"
#include "purge.h"
//#include "interdof.hpp"

//! \brief A object containing the list of forces
//! \author V. Richefeu
class Network
{
  struct particle_pair
  {
    body2d * i;
    body2d * j;
  };
  vector<inter2d*> linter_;
  vector<unsigned int> clist_;
 // vector<linterdof*>linterdof_;

  vector<particle_pair> superList_;
  vector<particle_pair> superListP_;

  vector<unsigned int> iNum_;
  vector<unsigned int> jNum_;

  double dsuperList_;
  double dsuperListP_;
  double dverlet_;
  bool useSuperList_;

  unsigned int speakLevel_;
	 
public:

  Network() 
    {
    dverlet_ = 1.0;
    speakLevel_ = 0;
    dsuperList_ = 1.0;
    dsuperListP_ = 1.0;
    useSuperList_ = false;
    }
  
  ~Network() { purge(linter_); }
  
 // vector<interdof*> & linterdof() {return linterdof_;}
  
  vector<inter2d*> &  linter() { return linter_; } 
  
  inter2d* inter(unsigned int k)   { return linter_[k]; }
  
  void read(istream &);
  //void read(istream &,Sample &);
  void write(ostream &);
  void write(const char *);
  void speak();
  
  void associate(Sample&);
  
  void clearTmpBodyNum()
    {
    iNum_.clear();
    jNum_.clear();
    }
  
  double & dverlet()       { return dverlet_; }
  double   dverlet() const { return dverlet_; }

  double & dsuperList()       { return dsuperList_; }
  double   dsuperList() const { return dsuperList_; }
  
  bool & useSuperList()       { return useSuperList_; }
  bool   useSuperList() const { return useSuperList_; }

  double & dsuperListP()       { return dsuperListP_; }
  double   dsuperListP() const { return dsuperListP_; }

vector<unsigned int>& clist() {return clist_;}

  unsigned int clist(unsigned int k) {return clist_[k];}
  //! \brief Build the list of potential interactions
  //! \author V. Richefeu
  void verlet(Sample*,  GroupRelationData*);
  
  //! \brief Build the super list...
  //! \author V. Richefeu
  void buildSuperList(Sample*,  GroupRelationData*);

  //! \brief Build the periodic super list...
  //! \author C. Voivret
  void buildSuperListP(Sample*,  GroupRelationData*);
  
  
  //! \brief Build the list of potential periodics interactions
  //! \author Charles Voivret
  void verletP(Sample*,  GroupRelationData*);
  
  //! Stock all the forces
  //! \author V. Richefeu
  void stock(vector<fsafe>&);
  
  //! Stock the forces belonging to 'sublist'
  //! \author V. Richefeu
  void stock(vector<fsafe>&, vector<unsigned int>& sublist);
  
  //! Retrieve the forces previously stocked
  //! \author V. Richefeu
  void retrieve(vector<fsafe>&);


};

#endif // _network_hpp
