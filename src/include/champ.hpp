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
		double resolution_ ; //resolution : definie par grid, ensemble de particules a prendre en compte autour d'un point
		double true_resolution_; // parametre de la ponderation expo, peut etre plus faible que resolution_ 
		string name_ ;
		typechamps type_;//identifie le champ pour le calcul
		double epsilon_;
	public:
		Champ();
		virtual ~Champ();
		virtual void calculMasse(const Grid&,Sample& ) {};
		virtual void writechamp(const Grid&) {};
		string getname() {return name_ ; }
		typechamps gettype() {return type_ ; }
		double ponderation(double);
};


class Champ_Scalaire : public Champ
{
	private:
	double * champ ;
	public :
	//nx,ny,nom,typechamp,resolution de la grille
	Champ_Scalaire(unsigned int,unsigned int, string,typechamps,double);
	~Champ_Scalaire();
	void calculMasse(const Grid&,Sample&);
	void writechamp(const Grid&);
};





#endif
