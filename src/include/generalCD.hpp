#ifndef _generalCD_h
#define _generalCD_h

#include <iostream>
#include <string>
#include "system.hpp"

class generalCD : public System
{

 protected:
	
 public:
  
  void read_parameters (istream&);
  void write_parameters(ostream&);
  void init ();
  void drive();
  void trans();
  void share();
  int  check();
  void stress_strain();

//void little_analyse(double );
  
  ~generalCD() { }
   generalCD(Sample* spl, Network* nwk, GroupRelationData * grpRel) : System(spl,nwk,grpRel) { }
   generalCD() : System() { }

};

#endif // _generalCD_h
