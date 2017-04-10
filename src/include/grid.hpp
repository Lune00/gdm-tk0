#ifndef _grid_hpp
#define _grid_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include"config.hpp"
#include"point_bis.hpp"
//Proprietes d'un point de la grille

class Grid{

	private:
		unsigned int nx_ , ny_ ;
		double xmin_, xmax_ , ymin_, ymax_ ;
		double dx_ , dy_ ;
		Point * array_ ; 
	public:
		Grid(){};
		~Grid();
		Grid(Config);
		double getX(int,int);
		double getY(int,int);
		Point * getarray() {return array_ ;}
		void writeGrid(string);
		void stockparticles(const Sample&);
};



Grid::Grid(Config parametres)
{
	nx_ = parametres.getnx() ;
	ny_ = parametres.getny() ;

	xmin_ = parametres.getxmin() ;
	xmax_ = parametres.getxmax() ;
	ymin_ = parametres.getymin() ;
	ymax_ = parametres.getymax() ;

	array_ = new Point [ nx_ * ny_ ];

	if(xmin_ > xmax_ || ymin_ > ymax_)
	{
		cerr<<"@Erreur definition de la metrique de la grille."<<endl;
		return ;
	}

	double Lx = xmax_ - xmin_ ;
	double Ly = ymax_ - ymin_ ;

	dx_ = Lx / nx_ ;
	dy_ = Ly / ny_ ;

	double x = xmin_;
	double y = ymin_;

	//Set coordinates:

	for(int j = 0 ; j!= ny_ ; j++)
	{
		for(int i = 0 ; i != nx_ ; i++)
		{
			array_[ i * ny_ + j ].setX(x);
			array_[ i * ny_ + j ].setY(y);
			x+=dx_;
		}
		y+=dy_;
		x = xmin_;
	}
	cerr<<"Metrics done."<<endl;
}


Grid::~Grid()
{
	delete[] array_;
}

double Grid::getX(int i,int j)
{
	return array_[ i * ny_ + j].getX();
}

double Grid::getY(int i,int j)
{
	return array_[ i * ny_ + j].getY();
}

void Grid::writeGrid(string filename)
{
	ofstream gridout (filename,ios::out);

	cout<<filename<<endl;
	for(int j = 0 ; j != ny_ ; j++)
	{
		for(int i = 0 ; i!= nx_ - 1 ; i++)
		{
			gridout<<getX(i,j)<<" "<<getY(i,j)<<" "<<dx_<<" 0."<<endl;
		}
	}

	for(int i = 0 ; i != nx_ ; i++)
	{
		for(int j = 0 ; j!= ny_ - 1 ; j++)
		{
			gridout<<getX(i,j)<<" "<<getY(i,j)<<" 0. "<<dy_<<endl; 
		}
	}

	gridout.close();
}	

//On stock les particules sur chaque point
void Grid::stockparticles(const Sample& spl)
{
	unsigned int N = spl.lbody().size();
	cerr<<"Nombre de particules a stocker : "<< N <<endl;



}
#endif //_grid_hpp
