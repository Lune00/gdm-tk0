#ifndef _cohesionLaw_h
#define _cohesionLaw_h

#include "inter2d.hpp"


class cohesionLaw
{

protected:

	//inter2d * inter;
	string name_;

public:

	virtual double  fco(  inter2d *)     = 0;  // cohesion force
	virtual double  dAct( inter2d *)     = 0;  // Activation distance
	virtual void    read(istream &)      = 0;
	virtual ~cohesionLaw() {};
	
	static cohesionLaw * factory(string type);
	
	string name() {return( name_);}

};

#endif // _cohesionLaw_h

