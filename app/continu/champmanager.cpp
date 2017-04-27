#include"champmanager.hpp"

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
	double resolution = parametres.getl();

	cerr<<"Initialisation des champs:"<<endl;

	if(parametres.getcalcmasse()) 
	{
		Champ * masse = new Champ_Scalaire(nx,ny,"masse",t_masse,resolution) ;
		lchamps_.push_back(masse);
		
	}
	if(parametres.getcalcmomentum()) 
	{
		Champ * momentum = new Champ_Vectoriel(nx,ny,"momentum",t_momentum,resolution) ;
		lchamps_.push_back(momentum);
	}
	if(parametres.getcalcvitesse() && parametres.getcalcmasse() && parametres.getcalcmomentum()) 
	{
		Champ * vitesse = new Champ_Vectoriel(nx,ny,"vitesse",t_vitesse,resolution) ;
		lchamps_.push_back(vitesse);
	}
}

void ChampManager::calculChamps(const Grid& grid,Sample& spl)  
{

	for (std::vector<Champ*>::iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		switch ((*it)->gettype())
		{
			case t_masse :
				(*it)->calculMasse(grid,spl);
				break;
			case t_momentum : 
				(*it)->calculMomentum(grid,spl);
				break;
			case t_vitesse : 
				(*it)->calculVitesse(getchamp("masse"),getchamp("momentum"));
				break;
		}
	}
}


const Champ* ChampManager::getchamp(string name) const
{
	for (std::vector<Champ*>::const_iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		if((*it)->getname() == name) return (*it);
	}
	return NULL;
}


void ChampManager::writeChamps(const Grid& grid)
{
	for (std::vector<Champ*>::iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		switch ((*it)->gettype())
		{
			case t_masse :
				(*it)->writechamp(grid);
				break;
			case t_momentum : 
				(*it)->writechamp(grid);
				break;
			case t_vitesse : 
				(*it)->writechamp(grid);
				break;
		}
	}
}
