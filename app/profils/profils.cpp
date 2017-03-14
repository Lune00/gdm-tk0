#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"








void profilMoyenVitesse(unsigned int Nbins)
{
	double vx[Nbins];
	double vxsquare[Nbins];
	double y[Nbins];

	if( Nbins == 0 ) {
		cerr<<"@profilMoyenVitesse le nombre de Nbins ne peut pas etre nul."<<endl;
		return;
	}

	ifstream is("Analyse/SProfile.txt");
	if(!is.is_open()) {
		cerr<<"@profilMoyenVitesse impossible ouvrir SProfile.txt"<<endl;
		return;
	}
	else
	{
		cout << k<< endl;
		double time, yt , vxt ;
		unsigned int j = 0 ;
		while(is)
		{
			if( j % Nbins == 0 ) j = 0 ;
			is >> time >> vxt >> yt ; 
			//cout<<time<<" "<<vxt<<endl;
			vx[j] += vxt ;
			vxsquare[j] += vxt * vxt ;
			y[j] = yt ;
			j++;
		}
		is.close();

		for (unsigned int i = 0 ; i < Nbins ; i ++ )
		{
			vx[i]       /= (double) Nbins;
			vxsquare[i] /= (double) Nbins;
			vxsquare[j] -= vx[i] * vx[i] ;
		}


		//Normalisation par vitesse de la plaque / par h
		//Sortir a la fois brut et normalise
		//y/h vx/v sqrt(dvx)/v y vx dvx

	}


}
int main (int argc,char **argv)
{


	bool calcProfilVitesse = false;
	unsigned int NbinsV = 0 ;

	string token;

	ifstream is(argv[1]);

	is>>token;
	while(is)
	{
		if(token=="vitesse"){
			calcProfilVitesse = true ;
			is>>NbinsV;
		}
		is>>token;
	}

	if(calcProfilVitesse) {
		cout<<".Calcul profil moyen de vitesse - Bins : "<<NbinsV<<endl;
		profilMoyenVitesse(NbinsV);
	}


	return 0;
}


