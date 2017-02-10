//! \class MDalgo
//! \brief Routines for Molecular Dynamics
//! \author Vincent Richefeu

#ifndef _MDalgo_h
#define _MDalgo_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "algo.hpp"
#include "near.hpp"
#include "MDforceLaw.hpp"


class MDalgo : public Algo
{	
  MDforceLaw * law_;

  double        dt_;      //!< Time increment
  double        tm_;      //!< Cumulative time
  unsigned int  nwpost_;  //!< Number of steps between 2 records
  unsigned int  nver_;    //!< Number of steps between 2 updates of the Verlet list
  unsigned int  nsi_;     //!< Initial step number
  unsigned int  nsf_;     //!< Final step number
  unsigned int  ns_;      //!< Current step number
		
public:

  MDalgo(Sample* spl, Network* nwk, System* sys, GroupData* grpDat, GroupRelationData* grpRel, MDforceLaw * law) 
    : Algo(spl, nwk, sys, grpDat, grpRel), law_(law)
    { tm_ = 0.0; nver_ = 1; } 
  
  MDalgo() { tm_ = 0.0; nver_ = 1; }
	
  void read_parameters  (istream&);
  void write_parameters (ostream&);

  //! Initialize the computation
  void stand();
  
  //! Update positions and velocities for one step
  void velocityVerletStep(unsigned int ns);
  
  //! ...
  void forces();
  
  //! ...
  void hand();

  unsigned int contactStatut( inter2d *);

  //! Return true if some post-treatments have to be done 
  bool wpost() { return (ns_%nwpost_ == 0); }

  unsigned int    nsi()     const { return nsi_; }
  unsigned int  & nsi()           { return nsi_; }
  unsigned int    ns()      const { return ns_; }
  unsigned int  & ns()            { return ns_; }
  unsigned int    nsf()     const { return nsf_; }
  unsigned int  & nsf()           { return nsf_; }
  double          dt()      const { return dt_; }
  double        & dt()            { return dt_; }
  unsigned int    nwpost()  const { return nwpost_; }
  unsigned int  & nwpost()        { return nwpost_; }
  unsigned int    nver()    const { return nver_; }
  unsigned int  & nver()          { return nver_; }
  double          tm()      const { return tm_; }
  double        & tm()            { return tm_; }
	
};

#endif // _MDalgo_h
