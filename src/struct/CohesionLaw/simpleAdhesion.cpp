#ifndef _simpleAdhesion_h
#define _simpleAdhesion_h

#include "cohesionLaw.hpp"
#define MIN(A,B) ((A)<(B) ? (A):(B))

class simpleAdhesion : public cohesionLaw
{
protected: 
	double force_;
public:
	double fco(  inter2d * inter ) {return force_;} 
	double dAct( inter2d * inter ) {return 0.;}
	void   read( istream & in )    { in >> force_;}
	simpleAdhesion():force_(0){ name_="simple";}
	
	~simpleAdhesion() {}
};

#endif // _weightedAdhesion_h




