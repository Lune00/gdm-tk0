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
		double resolution_ ;
		Point * array_ ; 
		std::vector<voisin> motif_ ;
	public:
		Grid(){};
		Grid(double,double,double,int); //init a grid with dx,dy,l_ (motif)
		~Grid();
		Grid(Config);
		double getX(int,int);
		double getY(int,int);
		double getdx() {return dx_ ; }
		double getdy() {return dy_ ; }
		Point * getarray() {return array_ ;}
		void writeGrid(string);
		void setcoordinates();
		void initmotif(Config&);
		void repartition(const Sample&);
		int check(Config);
		bool recouvrement(Point,Point);
		Point getPoint(int,int);
		bool motifinit();
};


#endif //_grid_hpp
