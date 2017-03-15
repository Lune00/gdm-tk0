#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"
#include<algorithm>








void profilMoyenVitesse(unsigned int Nbins)
{
	double vx[Nbins];
	double vxsquare[Nbins];
	double y[Nbins];
	std::fill_n(vx,Nbins,0.);
	std::fill_n(vxsquare,Nbins,0.);


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
		double time, yt , vxt ;
		unsigned int j = 0 ;
		double tr;
		while(is)
		{
			if( j % Nbins  == 0 ) j = 0 ;

			is >> time >> yt >> vxt >> tr ; 

			if(is.eof()) break;

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

		double v ;
		ifstream printSystem("Analyse/system.txt");

		if(!printSystem.is_open()){
			cerr<<"@profilMoyenVitesse fichier Analyse/system.txt manquant pour la normalisation"<<endl;
			return;
		}

		while(printSystem)
		{
			printSystem >> tr >> tr >> tr >> tr >> v >> tr >> tr ;
		}
		printSystem.close();
		cout<<"Vitesse de la plaque supérieure : "<<v<<endl;


		double ymin,ymax;
		ymin = y[0] ;
		ymax = y[Nbins-1];
		double h = ymax - ymin ;
		cout <<"h = "<<h<<endl;	
		//Output : 
		ofstream pv("profils/profvitesse.txt",ios::out);

		for (unsigned int i = 0 ; i < Nbins ; i++ )
		{
			pv << 2 * (y[i]-ymin) / h  - 1.<<" "<<vx[i]/v<<" "<<sqrt(vxsquare[i])/v<<" "<< y[i]<<" "<<vx[i]<<" "<<sqrt(vxsquare[i])<<endl;
		}
		pv.close();


	}


}
int main (int argc,char **argv)
{


	bool calcProfilVitesse = false;
	unsigned int NbinsV = 0 ;

	system("mkdir -p profils");
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


