#ifndef _algo_hpp
#define _algo_hpp

#include <iostream>
#include <fstream>
#include "sample.hpp"
#include "network.hpp"
#include "system.hpp"
#include "groupData.hpp"
#include "groupRelationData.hpp"

//! \brief Virtual class for algorithms
//! \author V. Richefeu
class Algo
{	

protected:

	Sample            * spl_;
	Network           * nwk_;
	System            * sys_;
	GroupData         * grpDat_;
	GroupRelationData * grpRel_;

	unsigned int speakLevel_;
	
	double dt_;  //!< Time increment (dynamics algorithm)

public:

	Algo(Sample* spl, Network* nwk, System* sys, GroupData* grpDat, GroupRelationData* grpRel) 
		: spl_(spl), nwk_(nwk), sys_(sys), grpDat_(grpDat), grpRel_(grpRel){ } 

	Algo() 
	{
		speakLevel_ = 1;
		spl_    = 0;
		nwk_    = 0;
		sys_    = 0;
		grpDat_ = 0;
		grpRel_ = 0;
	}

	virtual ~Algo() { }

	void plugSample   (Sample   * spl) { spl_ = spl; }
	void plugNetwork  (Network  * nwk) { nwk_ = nwk; }
	void plugSystem   (System   * sys) { sys_ = sys; }
	void plugGroupData (GroupData * grpDat) { grpDat_ = grpDat; }  
	void plugGroupRelationData (GroupRelationData * grpRel) { grpRel_ = grpRel; }
	
	void algoFill( );

	virtual void read_parameters  (istream&) = 0;
	virtual void write_parameters (ostream&) = 0;

	virtual void stand()            = 0;
	virtual void look()             = 0;
	virtual void hand(unsigned int) = 0;
	virtual void step()             = 0;
	virtual void speak()            = 0;
	
	virtual unsigned int contactStatut( inter2d *)=0;
	
	static Algo* factory(string type);
	
	double          dt()      const { return dt_; }
	double        & dt()            { return dt_; }
};

#endif // _algo_hpp

