#ifndef _weightedAdhesion_h
#define _weightedAdhesion_h

#include "cohesionLaw.hpp"


class weightedAdhesion : public cohesionLaw
{
protected: 
	double AdEnergy_;
public:
	double fco(  inter2d * inter ) {return AdEnergy_;} 
	double dAct( inter2d * inter ) {return 0.;}
	void   read( istream & in )    { in >> AdEnergy_;}
	weightedAdhesion() : AdEnergy_(0){ name_="weight";}
	
	~weightedAdhesion() {}
};

#endif // _weightedAdhesion_h




