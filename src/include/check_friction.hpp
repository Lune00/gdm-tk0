#ifndef _Check_friction_h
#define _Check_friction_h

#include <iostream>
#include <string>
#include "system.hpp"

#include "circularProbe.hpp"
#include "stress.hpp"
#include "tensor.hpp"

//! \brief A biaxial system. Left and Bottom walls are fixed and the
//!  other walls can be submitted to VELOCITIES, FORCES or PRESSURES. 
//! \author V. Richefeu
class Check_friction : public System
{

protected:

	unsigned int xmod_, ymod_;
	double       xval_, yval_;

	polyg * one_;
	polyg * two_;
	polyg * three_;

	
	bool gravity_;

public:

	void read_parameters (istream&);
	void write_parameters(ostream&);
	void init ();
	void drive();
	void trans();
	void share();
	int  check();
	
	polyg * one()  {return one_;}
	polyg * two()  {return two_;}
//	polyg * three()  {return left_;}
//	void little_analyse(double);
  	void stress_strain();


	~Check_friction() { }

	Check_friction(Sample* spl, Network* nwk, GroupRelationData * grpRel) : System(spl,nwk,grpRel)  
	{ 
		xval_ = yval_ = 0.0;
		xmod_ = ymod_ = _PRESSURE;
		gravity_ = false;
	}

	Check_friction() : System() 
	{ 
		xval_ = yval_ = 0.0;
		xmod_ = ymod_ = _VELOCITY;
		gravity_ = false;
	}

// TODO: move this function outside (in namespace gdm)
	unsigned int loadingMode(string mode)
	{
		if      (mode == "PRESSURE") return _PRESSURE;
		else if (mode == "VELOCITY") return _VELOCITY;
		else if (mode == "FORCE")    return _FORCE;
		else 
		{
			cout<<" loading Mode : bad parameter "<<endl;
		    return _VELOCITY;
		}
	}
	
// TODO: move this function outside (in namespace gdm)
	string loadingMode(unsigned int which)
	{
		string mode("_VELOCITY");
		if (which == _PRESSURE) mode = "PRESSURE";
		if (which == _VELOCITY) mode = "VELOCITY";
		if (which == _FORCE)    mode = "FORCE";
		return mode;
	}
};

#endif // _Check_friction_h
