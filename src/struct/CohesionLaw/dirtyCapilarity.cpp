#ifndef _weightedAdhesion_h
#define _weightedAdhesion_h

#include "cohesionLaw.hpp"


class dirtyCapilarity : public cohesionLaw
{
protected: 
	double force_;
	double debondingDistance_;
	
	
public:
	double fco(  inter2d * inter ) {return force_;} 
	double dAct( inter2d * inter ) {return debondingDistance_;}//peut etre ponderee par les diam des part en interaction...
	void   read( istream & in )    { in >> force_;}
	
	weightedAdhesion() : force_(0),debondingDistance_(0){ name_="dirtyCapilarity";}
	
	~weightedAdhesion() {}
};

#endif // _weightedAdhesion_h




