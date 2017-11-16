#include "config.hpp"

using namespace std;

Config::Config()
{
	nx = 0 ;
	ny = 0 ;
	istart = 0 ;
	iend = 0 ;
	di = 0 ;
	xmin = 0. ;
	xmax = 0. ;
	ymin = 0. ;
	ymax = 0. ;
	metrics = "";
	fichsim = "Simu.sim";
	calc_masse = false ;
	calc_momentum = false ;
	l_ = 0. ;
	dossierparent = "continu";
}

//Initialisation uniquement ici des parametres globaux de l'analyse
void Config::init(ifstream& is)
{
	if(!is)
	{
		cerr << "@initContinu : cannot open file " << endl;
		return;
	}

	string token ;
	is >> token ;
	while(is)
	{
		if(token=="Grid{")
		{
			is >> token;
			while(is)
			{
				if(token=="nx") is >>nx ; 
				if(token=="spatial") is >> l_ ; 
				if(token=="ny") is >>ny ; 
				if(token=="metrics") 
				{
					is >> token;
					readMetrics(token);
				}
				else if (token=="}") break;
				is >> token;
			}
		}
		if(token=="Data{")
		{
			is >> token;
			while(is)
			{
				if(token=="ini") is >>istart ; 
				if(token=="end") is >>iend ; 
				if(token=="di") is >>di; 
				if(token=="Simu") is >> fichsim;
				else if (token=="}") break;
				is >> token;
			}
		}
		if(token=="Champs{")
		{
			is >> token;
			while(is)
			{
				if (token=="masse") calc_masse = true ;
				if (token=="momentum") calc_momentum = true ;
				if (token=="vitesse") calc_vitesse = true ;
				if (token=="}") break;
				is >> token;
			}
		}

		is >> token;
	}

	//Initialisation de la résolution spatiale initiale:
	{
		l_ *= rmean_ ;
		cerr<<"---------------------------------"<<endl;
		cerr<<"initial spatial resolution (l_) = "<<l_<<endl;
		cerr<<"initial spatial resolution / average radius = "<<l_/rmean_<<endl;
		cerr<<"---------------------------------"<<endl;
	}

	string commande = "mkdir -p " + dossierparent ;
	std::system(commande.c_str());
}

void Config::readMetrics(string file)
{
	cerr<<"Fichier init metrics : "<<file<<endl;
	ifstream is(file.c_str());
	string token;
	is >> token;
	while(is)
	{
		if(token=="xmin") is >> xmin;
		if(token=="xmax") is >> xmax;
		if(token=="ymin") is >> ymin;
		if(token=="ymax") is >> ymax;
		if(token=="rmean") is >> rmean_ ;
		if(token=="bw") is >> bandwidth ;
		is >> token;
	}
}
