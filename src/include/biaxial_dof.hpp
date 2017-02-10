#ifndef _biaxial_dof_h
#define _biaxial_dof_h

#include <iostream>
#include <string>
#include "system.hpp"

#include "circularProbe.hpp"
#include "stress.hpp"
#include "tensor.hpp"

//! \brief A biaxial_dof system. Left and Bottom walls are fixed and the
//!  other walls can be submitted to VELOCITIES, FORCES or PRESSURES. 
//! \author V. Richefeu
class Biaxial_dof : public System
{

protected:

	unsigned int xmod_, ymod_;
	double       xval_, yval_;

	rline * bottom_;
	rline * top_;
	rline * left_;
	rline * right_;

	unsigned int idBottom_,idTop_,idLeft_,idRight_;
	bool adjust_;
	bool gravity_;

public:

	void read_parameters (istream&);
	void write_parameters(ostream&);
	void init ();
	void drive();
	void trans();
	void share();
	int  check();
	
	rline * bottom()  {return bottom_;}
	rline * top()  {return top_;}
	rline * left()  {return left_;}
	rline * right()  {return right_;}
   void stress_strain();

//	void little_analyse(double);

	~Biaxial_dof() { }

	Biaxial_dof(Sample* spl, Network* nwk, GroupRelationData * grpRel) : System(spl,nwk,grpRel)  
	{ 
		xval_ = yval_ = 0.0;
		xmod_ = ymod_ = _PRESSURE;
		idBottom_ = 0;
		idTop_    = 1;
		idLeft_   = 2;
		idRight_  = 3;
		adjust_   = false;
		gravity_ = false;
	}

	Biaxial_dof() : System() 
	{ 
		xval_ = yval_ = 0.0;
		xmod_ = ymod_ = _VELOCITY;
		idBottom_ = 0;
		idTop_    = 1;
		idLeft_   = 2;
		idRight_  = 3;
		adjust_   = false;
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

#endif // _biaxial_dof_h
