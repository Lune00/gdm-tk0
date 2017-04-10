#ifndef _config_hpp
#define _config_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>

using namespace std;

class Config{

	private:

		int nx ;
		int ny ;
		int istart, iend, di;
		double xmin, xmax, ymin, ymax ;
		//Resolution spatiale:
		double l_ ;
		string metrics;
		string fichsim;
		bool calc_masse;

	public:
		Config();
		~Config(){};
		void init(ifstream&);
		void readMetrics(string);

		int getnx() {return nx ;}
		int getny() {return ny ;}
		int getistart() {return istart;}
		int getiend() {return iend;}
		int getdi() {return di;}

		double getxmin() {return xmin ;}
		double getymin() {return ymin ;}
		double getxmax() {return xmax ;}
		double getymax() {return ymax ;}
		string getfichsim() {return fichsim;}

		bool getcalcmasse() {return calc_masse;}
		double getl() {return l_ ; }
};


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
				if (token=="}") break;
				is >> token;
			}
		}

		is >> token;
	}
}

void Config::readMetrics(string file)
{
	ifstream is(file);
	string token;
	is >> token;
	while(is)
	{
		if(token=="xmin") is >> xmin;
		if(token=="xmax") is >> xmax;
		if(token=="ymin") is >> ymin;
		if(token=="ymax") is >> ymax;
		is >> token;
	}
}




#endif
