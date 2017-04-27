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
#include "bitmap_image.hpp"

//Differents types de champs calculables:
enum typechamps{t_masse,t_momentum,t_vitesse};

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
		virtual void calculMomentum(const Grid&,Sample& ) {};
		virtual void calculVitesse(const Champ*,const Champ*) {};
		virtual void writechamp(const Grid&) {};
		//J'ai pas trouve mieux, ou alors formater le retour de valeur du champ
		//Par exemple une struct valeur avec 1,2,4 valeurs selon vec,tens,sca
		virtual double getchamp_ij(int,int) const   {return 0.; }
		virtual double getchamp_ij_x(int,int) const {return 0.; };
		virtual double getchamp_ij_y(int,int) const {return 0.; };
		string getname() {return name_ ; }
		typechamps gettype() {return type_ ; }
		double ponderation(double);
		virtual void drawchamp() {};
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
		void drawchamp();
		double getchamp_ij(int i,int j) const {return champ[i*ny_+j];}
};

class Champ_Vectoriel : public Champ
{
	private:
		double * champx;
		double * champy;
	public:
		Champ_Vectoriel(unsigned int, unsigned int, string, typechamps, double);
		~Champ_Vectoriel();
		void calculMomentum(const Grid&,Sample&);
		void calculVitesse(const Champ*,const Champ*) ;
		void writechamp(const Grid&);
		double getchamp_ij_x(int i,int j) const {return champx[i*ny_+j];}
		double getchamp_ij_y(int i,int j) const {return champy[i*ny_+j];}
};



#endif
