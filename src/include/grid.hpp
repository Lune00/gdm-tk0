#ifndef _grid_hpp
#define _grid_hpp

#include <iostream>
#include <set>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "config.hpp"
#include "point_bis.hpp"
#include "sample.hpp"

class Grid{

	private:
		unsigned int nx_ , ny_ ;
		unsigned int nb_ ;
		double xmin_, xmax_ , ymin_, ymax_ ;
		double bandwitdh_;
		double dx_ , dy_ ;
		double resolution_ ;
		Point * array_ ; 
		std::vector<Point> motif_ ;
	public:
		Grid(){};
		Grid(double,double,double,int); //init a grid with dx,dy,l_ (motif)
		~Grid();
		Grid(Config);
		double getX(int,int) const;
		double getY(int,int) const;
		double getnx() {return nx_ ;}
		double getny() {return ny_ ;}
		double getdx() {return dx_ ; }
		double getdy() {return dy_ ; }
		double getbandwidth() {return bandwitdh_;}
		double getResolution() {return resolution_;}
		Point * getarray() {return array_ ;}
		void writeGrid(string);
		void setcoordinates();
		void initmotif(Config&);
		void repartition( Sample&);
		int check(Config);
		bool recouvrement(Point,Point);
		Point getPoint(int,int);
		Point& returnPoint (int,int);
		Point& readPoint (int,int) const;
		bool motifinit();
		void updatePoints(int,int,body2d*);
		void clearPoints();
		double getdistance(Point,body2d*);
		bool belongtopoint(Point,body2d*);
		bool out(int,int);
};


#endif //_grid_hpp
