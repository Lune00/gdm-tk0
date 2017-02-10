#ifndef _nocohesion_h
#define _nocohesion_h

#include "cohesionLaw.hpp"


class nocohesion : public cohesionLaw
{
protected: 

public:
	double fco(  inter2d * inter) {return 0.;} 
	double dAct( inter2d * inter) {return 0.;}
	void   read( istream & in ){ }
	
	nocohesion(){ name_="nocohesion"; }
	
	~nocohesion() {}
};

#endif // _nocohesion_h




