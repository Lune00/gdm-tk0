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


void profilMoyenStress(unsigned int Nbins)
{
	double xy[Nbins];
	double yy[Nbins];
	double y[Nbins];
	std::fill_n(xy,Nbins,0.);
	std::fill_n(yy,Nbins,0.);


	if( Nbins == 0 ) {
		cerr<<"@profilMoyenVitesse le nombre de Nbins ne peut pas etre nul."<<endl;
		return;
	}

	ifstream is("Analyse/StressProfile.txt");
	if(!is.is_open()) {
		cerr<<"@profilMoyenStress impossible ouvrir StressProfile.txt"<<endl;
		return;
	}
	else
	{
		double time, yt , xy_s, yy_s ;
		unsigned int j = 0 ;
		unsigned int tic = 0 ;
		double tr;
		//?
		while(is)
		{
			if( j % Nbins  == 0 ){ j = 0 ; tic++ ; }

			is >> time >> yt >> tr >> xy_s >> yy_s >> tr >> tr ; 

			if(is.eof()) break;

			xy[j] +=  xy_s * (-1.)  ;
			yy[j] +=  yy_s  ;
			y[j] = yt ;
			j++;
		}

		tic--;
		is.close();

		for (unsigned int i = 0 ; i != Nbins ; i ++ )
		{
			xy[i]       /= (double) tic;
			yy[i]       /= (double) tic;
		}

		double ymin,ymax;
		ymin = y[0] ;
		ymax = y[Nbins-1];
		double h = ymax - ymin ;
		cout <<"h = "<<h<<endl;	
		//Output : 
		ofstream pv("profils/profstress.txt",ios::out);
		pv<<"# 2y/h-1 mu musquare"<<endl;

		for (unsigned int i = 0 ; i != Nbins ; i++ )
		{
			//	cout<<i<<" "<<vxsquare[i]<<" "<<vx[i]<<" "<<endl;
			pv << 2 * (y[i]-ymin) / h  - 1.<<" "<<xy[i]<<" "<<yy[i]<<endl;
		}
		pv.close();


	}
}


void profilMoyenShearRate(unsigned int Nbins)
{
	double vx[Nbins];
	double vxsquare[Nbins];
	double y[Nbins];
	std::fill_n(vx,Nbins,0.);
	std::fill_n(vxsquare,Nbins,0.);


	if( Nbins == 0 ) {
		cerr<<"@profilMoyenShearRate le nombre de Nbins ne peut pas etre nul."<<endl;
		return;
	}

	ifstream is("Analyse/ShearProfile.txt");
	if(!is.is_open()) {
		cerr<<"@profilMoyenShearRate impossible ouvrir ShearProfile.txt"<<endl;
		return;
	}
	else
	{
		double time, yt , vxt ;
		unsigned int j = 0 ;
		unsigned int tic = 0 ;
		double tr;
		//?
		while(is)
		{
			if( j % Nbins  == 0 ){ j = 0 ; tic++ ; }

			is >> time >> yt >> vxt  ; 

			if(is.eof()) break;

			vx[j] += vxt ;
			vxsquare[j] += vxt * vxt ;
			y[j] = yt ;
			j++;
			if(j==0) cout<<vx[0]<<" "<<vxsquare[0]<<endl;
		}

		tic--;
		is.close();

		for (unsigned int i = 0 ; i != Nbins ; i ++ )
		{
			vx[i]       /= (double) tic;
			vxsquare[i] /= (double) tic;
			vxsquare[i] -= vx[i] * vx[i] ;
		}


		//Normalisation par le taux de deformation applique
		//Sortir a la fois brut et normalise
		//y/h shearrate/shearRateImpose fluctuations/shearRateImpose

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
		ofstream pv("profils/profilshearrate.txt",ios::out);
		pv<<"# 2y/h-1 vx/v dvx/v y vx dvx"<<endl;

		//Si symetrique !!!!
		double shearRateImpose = 2 * v / h ;
		for (unsigned int i = 2 ; i != Nbins-2 ; i++ )
		{
			pv << 2 * (y[i]-ymin) / h  - 1.<<" "<<vx[i]/shearRateImpose<<" "<<sqrt(vxsquare[i])/shearRateImpose<<" "<< y[i]<<" "<<vx[i]<<endl;
		}
		pv.close();

	}


}

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
		unsigned int tic = 0 ;
		double tr;
		//?
		while(is)
		{
			if( j % Nbins  == 0 ){ j = 0 ; tic++ ; }

			is >> time >> yt >> vxt >> tr ; 

			if(is.eof()) break;

			vx[j] += vxt ;
			vxsquare[j] += vxt * vxt ;
			y[j] = yt ;
			j++;
			if(j==0) cout<<vx[0]<<" "<<vxsquare[0]<<endl;
		}

		tic--;
		is.close();

		for (unsigned int i = 0 ; i != Nbins ; i ++ )
		{
			vx[i]       /= (double) tic;
			vxsquare[i] /= (double) tic;
			vxsquare[i] -= vx[i] * vx[i] ;
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
		pv<<"# 2y/h-1 vx/v dvx/v y vx dvx"<<endl;

		for (unsigned int i = 2 ; i != Nbins-2 ; i++ )
		{
			//	cout<<i<<" "<<vxsquare[i]<<" "<<vx[i]<<" "<<endl;
			pv << 2 * (y[i]-ymin) / h  - 1.<<" "<<vx[i]/v<<" "<<sqrt(vxsquare[i])/v<<" "<< y[i]<<" "<<vx[i]<<" "<<sqrt(vxsquare[i]) / (double) tic<<endl;
		}
		pv.close();


	}


}

void profilMoyenTemperature(unsigned int Nbins)
{
	double vx[Nbins];
	double vxsquare[Nbins];
	double y[Nbins];
	std::fill_n(vx,Nbins,0.);
	std::fill_n(vxsquare,Nbins,0.);


	if( Nbins == 0 ) {
		cerr<<"@profilMoyenTemperature le nombre de Nbins ne peut pas etre nul."<<endl;
		return;
	}

	ifstream is("Analyse/TempProfile.txt");
	if(!is.is_open()) {
		cerr<<"@profilMoyenTemperature impossible ouvrir TempProfile.txt"<<endl;
		return;
	}
	else
	{
		double time, yt , vxt ;
		unsigned int j = 0 ;
		unsigned int tic = 0 ;
		double tr;
		//Pour l'instant on fait le profil de la trace
		while(is)
		{
			if( j % Nbins  == 0 ){ j = 0 ; tic++ ; }

			is >> time >> tr >> tr >> yt >> vxt >> tr >> tr >> tr >> tr ; 

			if(is.eof()) break;

			vx[j] += vxt ;
			vxsquare[j] += vxt * vxt ;
			y[j] = yt ;
			j++;
			if(j==0) cout<<vx[0]<<" "<<vxsquare[0]<<endl;
		}

		tic--;
		is.close();

		for (unsigned int i = 0 ; i != Nbins ; i ++ )
		{
			vx[i]       /= (double) tic;
			vxsquare[i] /= (double) tic;
			vxsquare[i] -= vx[i] * vx[i] ;
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
		double h = max(ymax,ymin) - min(ymax,ymin) ;
		cout <<"h = "<<h<<endl;	
		//Output : 
		ofstream pv("profils/proftemperature.txt",ios::out);
		pv<<"# 2y/h-1 trace dtrace"<<endl;

		for (unsigned int i = 2 ; i != Nbins-2 ; i++ )
		{
			//	cout<<i<<" "<<vxsquare[i]<<" "<<vx[i]<<" "<<endl;
			pv << 2 * (y[i]-ymin) / h  - 1.<<" "<<vx[i]<<" "<<sqrt(vxsquare[i]) / (double) tic<<endl;
		}
		pv.close();


	}
}

void profilMoyenZ(unsigned int Nbins)
{
	double vx[Nbins];
	double vxsquare[Nbins];
	double y[Nbins];
	std::fill_n(vx,Nbins,0.);
	std::fill_n(vxsquare,Nbins,0.);


	if( Nbins == 0 ) {
		cerr<<"@profilMoyenZ le nombre de Nbins ne peut pas etre nul."<<endl;
		return;
	}

	ifstream is("Analyse/ZProfile.txt");
	if(!is.is_open()) {
		cerr<<"@profilMoyenZ impossible ouvrir ZProfile.txt"<<endl;
		return;
	}
	else
	{
		double time, yt , vxt ;
		unsigned int j = 0 ;
		unsigned int tic = 0 ;
		//Pour l'instant on fait le profil de la trace
		while(is)
		{
			if( j % Nbins  == 0 ){ j = 0 ; tic++ ; }

			is >> time >> yt >> vxt  ; 

			if(is.eof()) break;

			vx[j] += vxt ;
			vxsquare[j] += vxt * vxt ;
			y[j] = yt ;
			j++;
			if(j==0) cout<<vx[0]<<" "<<vxsquare[0]<<endl;
		}

		tic--;
		is.close();

		for (unsigned int i = 0 ; i != Nbins ; i ++ )
		{
			vx[i]       /= (double) tic;
			vxsquare[i] /= (double) tic;
			vxsquare[i] -= vx[i] * vx[i] ;
		}


		//Normalisation par vitesse de la plaque / par h
		//Sortir a la fois brut et normalise
		//y/h vx/v sqrt(dvx)/v y vx dvx

		double ymin,ymax;
		ymin = y[0] ;
		ymax = y[Nbins-1];
		double h = ymax - ymin ;
		cout <<"h = "<<h<<endl;	
		//Output : 
		ofstream pv("profils/profZ.txt",ios::out);
		pv<<"# 2y/h-1 Z dZ"<<endl;

		for (unsigned int i = 2 ; i != Nbins-2 ; i++ )
		{
			//	cout<<i<<" "<<vxsquare[i]<<" "<<vx[i]<<" "<<endl;
			pv << 2 * (y[i]-ymin) / h  - 1.<<" "<<vx[i]<<" "<<sqrt(vxsquare[i]) <<endl;
		}
		pv.close();


	}
}

int main (int argc,char **argv)
{

	bool calcProfilVitesse = false;
	bool calcProfilStress = false;
	bool calcProfileShear = false;
	bool calcProfileTemp = false;
	bool calcProfileZ = false;
	unsigned int NbinsV = 0 ;
	unsigned int NbinsS = 0 ;
	unsigned int NbinsSh = 0 ;
	unsigned int NbinsST = 0 ;
	unsigned int NbinsZ = 0 ;

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
		if(token=="stress"){
			calcProfilStress = true ;
			is>>NbinsS;
		}
		if(token=="shearrate"){
			calcProfileShear = true ;
			is>>NbinsSh;
		}
		if(token=="temperature"){
			calcProfileTemp = true ;
			is>>NbinsST;
		}
		if(token=="Z"){
			calcProfileZ = true ;
			is>>NbinsZ;
		}
		is>>token;
	}

	if(calcProfilVitesse) {
		cout<<".Calcul profil moyen de vitesse - Bins : "<<NbinsV<<endl;
		profilMoyenVitesse(NbinsV);
	}
	if(calcProfilStress) {
		cout<<".Calcul profil moyen de stress (mu) - Bins : "<<NbinsS<<endl;
		profilMoyenStress(NbinsS);
	}
	if(calcProfileShear) {
		cout<<".Calcul profil moyen de shearrate - Bins : "<<NbinsSh<<endl;
		profilMoyenShearRate(NbinsSh);
	}
	if(calcProfileTemp) {
		cout<<".Calcul profil moyen de temperature - Bins : "<<NbinsST<<endl;
		profilMoyenTemperature(NbinsST);
	}
	if(calcProfileZ) {
		cout<<".Calcul profil moyen de Z - Bins : "<<NbinsZ<<endl;
		profilMoyenZ(NbinsZ);
	}


	return 0;
}


