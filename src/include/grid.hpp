#ifndef _grid_hpp
#define _grid_hpp
#include"config.hpp"
#include"champ.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
//Proprietes d'un point de la grille
class Point{

	private :
		double x_ , y_ ;
	public:
		Point(){x_ = 0. ; y_ = 0. ;}
		~Point(){};
		double getX(){return x_;};
		double getY(){return y_;};
		void setX(double x){x_ = x ;};
		void setY(double y){y_ = y ;};
};

class Grid{

	private:
		unsigned int nx_ , ny_ ;
		double xmin_, xmax_ , ymin_, ymax_ ;
		double dx_ , dy_ ;
		Point * array ; 
		std::vector<Champ*> lchamps_ ;
	public:
		Grid(){};
		~Grid();
		Grid(Config);
		void initChamps(Config);
		void calculChamps();
		double getX(int,int);
		double getY(int,int);
		void writeGrid(string);
};


#endif //_grid_hpp
Grid::Grid(Config parametres)
{
	nx_ = parametres.getnx() ;
	ny_ = parametres.getny() ;

	xmin_ = parametres.getxmin() ;
	xmax_ = parametres.getxmax() ;
	ymin_ = parametres.getymin() ;
	ymax_ = parametres.getymax() ;

	array = new Point [ nx_ * ny_ ];

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
			array[ i * ny_ + j ].setX(x);
			array[ i * ny_ + j ].setY(y);
			x+=dx_;
		}
		y+=dy_;
		x = xmin_;
	}
	cerr<<"Metrics done."<<endl;
	initChamps(parametres);
}

void Grid::initChamps(Config parametres)
{
	int nx = parametres.getnx();
	int ny = parametres.getny();
	cerr<<"Initialisation des champs:"<<endl;
	if(parametres.getcalcmasse()) 
	{
		Champ * masse = new Champ_Scalaire(nx,ny,"masse",t_masse) ;
		lchamps_.push_back(masse);
	}
}

Grid::~Grid()
{
	for (std::vector<Champ*>::iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		delete (*it);
	}
	delete[] array;
}

double Grid::getX(int i,int j)
{
	return array[ i * ny_ + j].getX();
}

double Grid::getY(int i,int j)
{
	return array[ i * ny_ + j].getY();
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

void Grid::calculChamps()
{

	for (std::vector<Champ*>::iterator it = lchamps_.begin(); it != lchamps_.end();it++)
	{
		switch ((*it)->gettype())
		{
			case t_masse : (*it)->calculMasse();
			case t_momentum : cerr<<"Coming."<<endl;
		}
	}
}
