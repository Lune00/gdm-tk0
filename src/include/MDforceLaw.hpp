//! \class MDforceLaw
//! \brief Force laws for Molecular Dynamics
//! \author Vincent Richefeu

#ifndef _MDforceLaw_h
#define _MDforceLaw_h

#include "network.hpp"
#include "GroupRelationData.hpp"

class MDforceLaw
{	
  
protected:
  
  Network  * nwk_;
  GroupRelationData * grpRel_;
		
public:
  
  MDforceLaw(Network* nwk, GroupRelationData* grpRel) 
    : nwk_(nwk), grpRel_(grpRel) { } 
  
  MDforceLaw() { }
	virtual ~MDforceLaw() { }
  
  void plugNetwork  (Network* nwk ) { nwk_ = nwk; }
  void plugGroupRelationData (GroupRelationData* grpRel) { grpRel_ = grpRel; }
  
  virtual void checkRequiredParameters() = 0;
  
  //virtual double fn   (unsigned int k) = 0;
  //virtual double ft   (unsigned int k) = 0;
  //virtual double frot (unsigned int k) = 0;
  
  virtual void computeForces(unsigned int k) = 0;
	
};

#endif // _MDforceLaw_h
