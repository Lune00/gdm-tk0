#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"








void profilMoyenVitesse(unsigned int Nbins,unsigned int Nstart,unsigned int Nend)
{
if( Nbins == 0 || (Nstart>Nend) ) {
	cerr<<"@profilMoyenVitesse erreur parametres."<<endl;
	return;
}



}
int main (int argc,char **argv)
{


	bool calcProfilVitesse = false;
	unsigned int NbinsV = 0 ;
	unsigned int Nstart = 0 ;
	unsigned int Nend = 0 ;

	string token;

	ifstream is(argv[1]);

	is>>token;
	while(is)
	{
		if(token=="vitesse"){
			calcProfilVitesse = true ;
			is>>NbinsV;
			is>>Nstart;
			is>>Nend;
		}
		is>>token;
	}

	if(calcProfilVitesse) {
		cout<<".Calcul profil moyen de vitesse - Bins : "<<NbinsV<<endl;
		profilMoyenVitesse(NbinsV,Nstart,Nend);
	}


	return 0;
}


