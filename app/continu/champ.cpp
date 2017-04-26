#include"champ.hpp"

using namespace std;

Champ::Champ(){}


Champ::~Champ(){}

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

void Champ_Scalaire::calculMasse(const Grid& grid,Sample& spl)
{

	cerr<<"On calcule."<<endl;
}

