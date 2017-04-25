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

//struct voisin{
//	int i_ ;
//	int j_;
//	voisin() : i_(0), j_(0) {} ;
//	voisin(int i,int j) : i_(i) , j_(j) {} ;
////	bool operator== (const voisin&b)
////	{
////		return ( i_ == b.i_ && j_ == b.j_);
////	}
//};
//
//
//struct comparevoisin{
//	bool operator() (const voisin&a,const voisin&b)
//	{
//		return ( ! ( a.i_ != b.i_ || a.j_ != b.j_ )) ;
//		//return ( ! (a.i_ == b.i_ && a.j_ == b.j_) ) ;
//	}
//};
class Grid{

	private:
		unsigned int nx_ , ny_ ;
		double xmin_, xmax_ , ymin_, ymax_ ;
		double dx_ , dy_ ;
		double resolution_ ;
		Point * array_ ; 
		std::vector<Point> motif_ ;
	public:
		Grid(){};
		Grid(double,double,double,int); //init a grid with dx,dy,l_ (motif)
		~Grid();
		Grid(Config);
		double getX(int,int);
		double getY(int,int);
		double getdx() {return dx_ ; }
		double getdy() {return dy_ ; }
		double getResolution() {return resolution_;}
		Point * getarray() {return array_ ;}
		void writeGrid(string);
		void setcoordinates();
		void initmotif(Config&);
		void repartition( Sample&);
		int check(Config);
		bool recouvrement(Point,Point);
		Point getPoint(int,int);
		Point& returnPoint(int,int);
		bool motifinit();
		void updatePoints(Point&,body2d*);
		double getdistance(Point,body2d*);
		bool belongtopoint(Point,body2d*);
		bool out(int,int);
};


#endif //_grid_hpp
