#ifndef _shear_h
#define _shear_h

#include <iostream>
#include <string>
#include "system.hpp"

//! \brief A simple shear system. 
//! \author V. Richefeu
class Shear : public System
{

protected:
	
	double pressure_;
	double rate_;
	
	double xc0_,yc0_;

	rline * bottom_;
	rline * top_;
	rline * left_;
	rline * right_;

	unsigned int idBottom_,idTop_,idLeft_,idRight_;
	bool adjust_;

public:

	void read_parameters (istream&);
	void write_parameters(ostream&);
	void init ();
	void drive();
	void trans();
	void share();
	int  check();
	//void little_analyse(double);
  	void stress_strain();

	

	~Shear() { }

	Shear(Sample* spl, Network* nwk, GroupRelationData * grpRel) : System(spl,nwk,grpRel) 
	{ 
		pressure_ = 1.0;
		rate_     = 0.0;
		idBottom_ = 0;
		idTop_    = 1;
		idLeft_   = 2;
		idRight_  = 3;
		adjust_   = false;
	}

	Shear() : System() 
	{ 
		pressure_ = 1.0;
		rate_     = 0.0;
		idBottom_ = 0;
		idTop_    = 1;
		idLeft_   = 2;
		idRight_  = 3;
		adjust_   = false;
	}

// TODO: move this function outside (in namespace gdm)
	unsigned int loadingMode(string mode)
	{
		if (mode == "PRESSURE") return _PRESSURE;
		if (mode == "VELOCITY") return _VELOCITY;
		if (mode == "FORCE")    return _FORCE;
		return _VELOCITY;
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

#endif // _shear_h
