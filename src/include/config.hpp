#ifndef _config_hpp
#define _config_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>


class Config{

	private:

		int nx ;
		int ny ;
		int istart, iend, di;
		double xmin, xmax, ymin, ymax ;
		//Resolution spatiale:
		double l_ ;
		double rmean_ ; //rayon moyen d'une particule de la sim
		std::string metrics;
		std::string fichsim;
		bool calc_masse;

	public:
		Config();
		~Config(){};
		void init(std::ifstream&);
		void readMetrics(std::string);

		int getnx() {return nx ;}
		int getny() {return ny ;}
		int getistart() {return istart;}
		int getiend() {return iend;}
		int getdi() {return di;}

		double getxmin() {return xmin ;}
		double getymin() {return ymin ;}
		double getxmax() {return xmax ;}
		double getymax() {return ymax ;}
		std::string getfichsim() {return fichsim;}

		bool getcalcmasse() {return calc_masse;}
		double getl() {return l_ ; }
};

#endif
