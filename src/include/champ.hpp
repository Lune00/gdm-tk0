#ifndef _champ_hpp
#define _champ_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "config.hpp"
#include "grid.hpp"

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
		virtual void calculMasse(const Grid&,Sample& ) {};
		string getname() {return name_ ; }
		typechamps gettype() {return type_ ; }
};


class Champ_Scalaire : public Champ
{
	private:
	double * champ ;
	public :
	Champ_Scalaire(unsigned int,unsigned int, string,typechamps);
	~Champ_Scalaire();
	void calculMasse(const Grid&,Sample&);
};





#endif
