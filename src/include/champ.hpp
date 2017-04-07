#ifndef _champ_hpp
#define _champ_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
// #include "grid.hpp"
#include "point_bis.hpp"

//Differents types de champs calculables:
enum typechamps{t_masse,t_momentum};


class Champ
{

	protected:
		unsigned int nx_ ;
		unsigned int ny_ ;
		string name_ ;
		typechamps type_;//identifie le champ pour le calcul
	public:
		Champ();
		virtual ~Champ();
		virtual void calculMasse(Simulation*,Point*) {};
		string getname() {return name_ ; }
		typechamps gettype() {return type_ ; }
};

Champ::Champ(){}


Champ::~Champ(){}


class Champ_Scalaire : public Champ
{
	private:
	double * champ ;
	public :
	Champ_Scalaire(unsigned int,unsigned int, string,typechamps);
	~Champ_Scalaire();
	void calculMasse(Simulation*,Point*);

};


Champ_Scalaire::Champ_Scalaire(unsigned int nx, unsigned int ny, string name, typechamps t)
{
	nx_ = nx ;
	ny_ = ny ;
	name_ = name ;
	champ = new double [ nx_ * ny_ ];
	type_ = t ;
}

Champ_Scalaire::~Champ_Scalaire()
{
	cerr<<"Champ scalaire "<<this->name_<<" libéré."<<endl;
	delete [] champ ;
}

void Champ_Scalaire::calculMasse(Simulation* mySimu, Point* grid)
{

	cerr<<"On calcule."<<endl;
}



#endif
