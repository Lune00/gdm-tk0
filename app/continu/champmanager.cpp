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
	cerr<<"Initialisation des champs:"<<endl;
	if(parametres.getcalcmasse()) 
	{
		Champ * masse = new Champ_Scalaire(nx,ny,"masse",t_masse) ;
		lchamps_.push_back(masse);
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
				cerr<<"Coming."<<endl;
				break;
		}
	}
}
