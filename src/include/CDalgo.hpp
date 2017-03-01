#ifndef _CDalgo_hpp
#define _CDalgo_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include "talk.hpp"
#include "algo.hpp"
#include "near.hpp"
#include "fsafe.hpp"
#include "Fragmentation.hpp"
#include "io.hpp" // pour network_write dans iter_2

//! \brief Routines for 2D Contact Dynamics
//! \author Vincent Richefeu
class CDalgo : public Algo
{	

	vector<fsafe> forcesSafe_;
	//vector<unsigned int> clist_; //!< List of contacts

	unsigned int  nitermn_; //!< Minimum number of iterations before check convergence
	unsigned int  niterconv_; //!< Number of converged iterations before break out
	unsigned int  nitermx_; //!< Maximum number of iterations
	unsigned int  niter_;   //!< Current number of iterations 
	unsigned int  nver_;    //!< Number of iterations between two update of the neighbor list
	unsigned int  nsuper_;    //!< Number of iterations between two update of the super neighbor list
	double        epsf_;    //!< Relative precision with respect to average normal force

	public:


	CDalgo(Sample* spl, Network* nwk, System* sys, GroupData* grpDat, GroupRelationData* grpRel) 
		: Algo(spl,nwk,sys,grpDat,grpRel)
	{ 
		nver_ = 1;
		niterconv_ = 1;
		nitermn_=1; 
		iterFunc = &CDalgo::iter_0; 
		speakLevel_ = 1;
	}

	CDalgo() 
	{
		nver_ = 1;
		nitermn_=1;
		niterconv_ = 1; 
		iterFunc = &CDalgo::iter_0; 
		speakLevel_ = 1;
	}

	~CDalgo() { }


	void read_parameters  (istream&);
	void write_parameters (ostream&);

	void speak();

	//! Initialize the computation
	void stand();

	//! Take care of contacts and forces
	void look();

	//! Increment time, Update Verlet list
	void hand(unsigned int);

	//! Update the contact list, local frames and relative velocities
	void contact();

	//! Compute force resultants
	void fres();

	//! Update positions and velocities for one step
	void step();

	//! Vérivier la condition résistance des contacts
	void resistance();

	//
	unsigned int contactStatut( inter2d *);

	unsigned int (CDalgo::*iterFunc) ();

	//! Gauss-Seidel Iteractions with contact, friction, and rolling friction
	unsigned int iter_0 (); 

	//! Gauss-Seidel Iteractions with Moreau's criteria (sup)
	unsigned int iter_1 (); 

	// Une etude speciale
	unsigned int iter_2 (); 

	unsigned int    nitermn() const { return nitermn_; }
	unsigned int  & nitermn()       { return nitermn_; }
	unsigned int    nitermx() const { return nitermx_; }
	unsigned int  & nitermx()       { return nitermx_; }
	unsigned int    niter()   const { return niter_; }
	unsigned int  & niter()         { return niter_; }
	unsigned int    nver()    const { return nver_; }
	unsigned int  & nver()          { return nver_; }
	double          epsf()    const { return epsf_; }
	double        & epsf()          { return epsf_; }
};

#endif // _CDalgo_hpp
