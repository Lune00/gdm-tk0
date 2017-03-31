#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"
#include<algorithm>


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
	public:
		Grid(){};
		~Grid();
		Grid(unsigned int, unsigned int);
		double getX(int,int);
		double getY(int,int);
		void initMetrics(double,double,double,double);
		void writeGrid(string);
};

//Ajouter un fichier a lire poux les min,max et set tous les paramètres
Grid::Grid(unsigned int nx , unsigned int ny)
{
	nx_ = nx ;
	ny_ = ny ;

	array = new Point [ nx_ * ny_ ];

	xmin_ = 0. ;
	ymin_ = 0. ;
	xmax_ = 0. ;
	ymax_ = 0. ;
	dx_ = 0. ;
	dy_ = 0. ;
}
Grid::~Grid()
{
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
void Grid::initMetrics(double xmin, double xmax, double ymin, double ymax)
{
	xmin_ = xmin ;
	xmax_ = xmax ;
	ymin_ = ymin ;
	ymax_ = ymax ;

	if(xmin_ > xmax_ || ymin_ > ymax_)
	{
		cerr<<"@Erreur definition de la metrique de la grille."<<endl;
		return ;
	}

	double Lx = xmax_ - xmin_ ;
	double Ly = ymax_ - ymin_ ;

	dx_ = Lx / nx_ ;
	dy_ = Ly / ny_ ;

	cerr<<"dx_ = "<<dx_<<endl ;
	cerr<<"dy_ = "<<dy_<<endl;

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

int main (int argc,char **argv)
{
	cout<<"Hello (continu) world !"<<endl;
	cout<<"Let's start."<<endl;
	Grid grid(128,128);
	grid.initMetrics(0.,10.,0.,10.);
	grid.writeGrid("grid.txt");
	return 0;
}


