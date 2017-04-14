#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"
#include "vecteur.hpp"
#include "grid.hpp"
#include "champ.hpp"


int main (int argc,char **argv)
{
	//Config utilisateur;
	//Gestionnaire de Champs;
	Config parametres;
	ChampManager MesChamps;

	ifstream is(argv[1]);
	//Initialisation des parametres:
	parametres.init(is) ;

	//Initialisation grille:
	Grid grid(parametres);
	//Grid tmp(grid.getdx(),grid.getdy(),parametres.getl());
	//grid.initmotif(tmp);
	//Check recouvrement:
	if(grid.check(parametres) == 1 ) return 0 ;

	grid.writeGrid("grid.txt");

	//Initialisation des champs:
	MesChamps.initChamps(parametres);

	return 0 ;

	Simulation * mySimu = new Simulation();
	mySimu->read_data(parametres.getfichsim().c_str());

	char nomFichier[100];

	for (unsigned int i = parametres.getistart() ; i != parametres.getiend(); i+=parametres.getdi()){

		sprintf(nomFichier,"spl_nwk/spl_nwk_%.4d.his",i);
		cerr<<"****** Chargement : "<<nomFichier<<endl;

		mySimu->load_history(nomFichier);
		mySimu->algo()->algoFill();

		//Calcul des champs
		MesChamps.calculChamps(mySimu,grid);
	}

	delete mySimu ;

	return 0;
}


