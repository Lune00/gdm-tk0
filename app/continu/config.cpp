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
	l_ = 0. ;
}

//Initialisation uniquement ici des parametres globaux de l'analyse
void Config::init(ifstream& is)
{
	cerr<<"))))))))))))))))))))))))))))))))))"<<endl;
	if(!is)
	{
		cerr << "@initContinu : cannot open file " << endl;
		return;
	}

	string token ;
	is >> token ;
	while(is)
	{
		cerr<<token<<"*"<<endl;
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
				if (token=="}") break;
				is >> token;
			}
		}

		is >> token;
	}
}

void Config::readMetrics(string file)
{
	cerr<<"Fichier init metrics : "<<file<<endl;
	ifstream is(file);
	string token;
	is >> token;
	while(is)
	{
		if(token=="xmin") is >> xmin;
		if(token=="xmax") is >> xmax;
		if(token=="ymin") is >> ymin;
		if(token=="ymax") is >> ymax;
		if(token=="rmean")
		{
			cerr<<"---------------------------------"<<endl;
			double mult ;
			is >> mult ;
			l_ *= mult ;
			cerr<<"l_ = "<<l_<<endl;
		}
		is >> token;
	}
}
