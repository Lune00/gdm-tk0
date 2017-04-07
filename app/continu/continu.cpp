#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"
#include <algorithm>
#include "vecteur.hpp"
#include "grid.hpp"


int main (int argc,char **argv)
{
	Config parametres;
	ifstream is(argv[1]);
	parametres.init(is) ;

	Grid grid(parametres);
	grid.writeGrid("grid.txt");

	Simulation * mySimu = new Simulation();
	mySimu->read_data(parametres.getfichsim().c_str());

	char nomFichier[100];

	for (unsigned int i = parametres.getistart() ; i != parametres.getiend(); i+=parametres.getdi()){

		sprintf(nomFichier,"spl_nwk/spl_nwk_%.4d.his",i);
		cerr<<"****** Chargement : "<<nomFichier<<endl;

		mySimu->load_history(nomFichier);
		mySimu->algo()->algoFill();

		//Calcul des champs
		grid.calculChamps(mySimu);
	}

	delete mySimu ;

	return 0;
}


