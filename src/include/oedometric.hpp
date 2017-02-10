#ifndef _oedometric_h
#define _oedometric_h

#include <iostream>
#include <string>
#include "system.hpp"

#include "circularProbe.hpp"
#include "rectangularProbe.hpp"
#include "stress.hpp"
#include "tensor.hpp"

//! \brief A biaxial system. Left and Bottom walls are fixed and the
//!  other walls can be submitted to VELOCITIES, FORCES or PRESSURES. 
//! \author V. Richefeu
class Oedometric : public System
{

protected:

	unsigned int k_; //time step counting
	double lx_;
	double xmin_,xmax_;
	double yval_;
	unsigned int ymod_;
	double ly_;
	double ymin_,ymax_;
	double y0_; // initial position of top wall
	double defy_; // axial strain
	
	double incr_; // percentage of charge increment

	rline * bottom_;
	rline * top_;

	unsigned int idBottom_,idTop_;
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
  	void stress_strain();

	
	rline * bottom()  {return bottom_;}
	rline * top()  {return top_;}
	
	
	

	~Oedometric() { }

	Oedometric(Sample* spl, Network* nwk) : System(spl,nwk) 
	{ 
		yval_ = 0.0;
		ymod_ = _VELOCITY;
		idBottom_ = 0;
		idTop_    = 1;
		adjust_   = false;
		gravity_ = false;

	}

	Oedometric() : System() 
	{ 
		yval_ = 0.0;
		ymod_ = _VELOCITY;
		idBottom_ = 0;
		idTop_    = 1;
		adjust_   = false;
		gravity_ = false;
		lx_ = 0.;
		xmin_ = xmax_ = 0.;
		ymin_ = ymax_ = 0.;
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

#endif // _biaxial_h
