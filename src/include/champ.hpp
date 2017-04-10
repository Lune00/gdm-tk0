#ifndef _champ_hpp
#define _champ_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "point_bis.hpp"
#include "config.hpp"

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


class ChampManager
{

	private:
		std::vector<Champ*> lchamps_ ;
	public:
		ChampManager();
		~ChampManager();
		void initChamps(Config);
		void calculChamps(Simulation*,Point*);
};

ChampManager::ChampManager(){}

ChampManager::~ChampManager()
{

	for (std::vector<Champ*>::iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		delete (*it);
	}
}

void ChampManager::initChamps(Config parametres)
{
	int nx = parametres.getnx();
	int ny = parametres.getny();
	cerr<<"Initialisation des champs:"<<endl;
	if(parametres.getcalcmasse()) 
	{
		Champ * masse = new Champ_Scalaire(nx,ny,"masse",t_masse) ;
		lchamps_.push_back(masse);
	}
}

void ChampManager::calculChamps(Simulation* mySimu, Point* array)
{

	for (std::vector<Champ*>::iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		switch ((*it)->gettype())
		{
			case t_masse :
				(*it)->calculMasse(mySimu,array);
				break;
			case t_momentum : 
				cerr<<"Coming."<<endl;
				break;
		}
	}
}

#endif
