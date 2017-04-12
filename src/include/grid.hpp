#ifndef _grid_hpp
#define _grid_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "config.hpp"
#include "point_bis.hpp"
#include "sample.hpp"


struct voisin{
	int i_ ;
	int j_;
};


class Grid{

	private:
		unsigned int nx_ , ny_ ;
		double xmin_, xmax_ , ymin_, ymax_ ;
		double dx_ , dy_ ;
		Point * array_ ; 
		std::vector<voisin> motif_ ;
	public:
		Grid(){};
		~Grid();
		Grid(Config);
		double getX(int,int);
		double getY(int,int);
		Point * getarray() {return array_ ;}
		void writeGrid(string);
		void initmotif(Config);
		void stockparticles(const Sample&);
};


#endif //_grid_hpp
