#ifndef _MDlaw0_hpp
#define _MDlaw0_hpp

#include <string>
#include "MDforceLaw.hpp"

//! \brief Elastique-linear contact and viscuous Coulomb friction
//! \author Vincent Richefeu
class MDlaw0 : public MDforceLaw
{	
		
public:
  
  MDlaw0(Network* nwk, GroupRelationData* grpRel) : MDforceLaw(nwk,grpRel) { } 
  MDlaw0() { }
  ~MDlaw0() { }
  
  void checkRequiredParameters();
  void computeForces(unsigned int k);
  /*
  inline double fn   (unsigned int k);
  inline double ft   (unsigned int k);
  inline double frot (unsigned int k);
	*/
};

#endif // _MDlaw0_hpp
