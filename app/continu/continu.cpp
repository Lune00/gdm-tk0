#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"
#include <algorithm>

struct Globals{
	int nx ;
	int ny ;
	int istart, iend, di;
	double xmin, xmax, ymin, ymax ;
	string metrics;
};

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
		Grid(Globals);
		double getX(int,int);
		double getY(int,int);
		void writeGrid(string);
};

//Ajouter un fichier a lire poux les min,max et set tous les paramètres
Grid::Grid(Globals parametres)
{
	nx_ = parametres.nx ;
	ny_ = parametres.ny ;

	xmin_ = parametres.xmin ;
	xmax_ = parametres.xmax ;
	ymin_ = parametres.ymin ;
	ymax_ = parametres.ymax ;

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


void readMetrics(string file, Globals& parametres)
{
	ifstream is(file);
	string token;
	is >> token;
	while(is)
	{
		if(token=="xmin") is >> parametres.xmin;
		if(token=="xmax") is >> parametres.xmax;
		if(token=="ymin") is >> parametres.ymin;
		if(token=="ymax") is >> parametres.ymax;
		is >> token;
	}
}
		
void initContinu(ifstream& is,Globals& parametres)
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
				if(token=="nx") is >>parametres.nx ; 
				if(token=="ny") is >>parametres.ny ; 
				if(token=="metrics") 
				{
					is >> token;
					readMetrics(token,parametres);
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
				if(token=="ini") is >>parametres.istart ; 
				if(token=="end") is >>parametres.iend ; 
				if(token=="di") is >>parametres.di; 
				else if (token=="}") break;
				is >> token;
			}
		}

		is >> token;
	}
}

int main (int argc,char **argv)
{
	Globals parametres;
	ifstream is(argv[1]);
	initContinu(is,parametres);
	Grid grid(parametres);
	cerr<<"nx = "<<parametres.nx<<" ny = "<<parametres.ny<<endl;
	cerr<<"Init metrics file : "<<parametres.metrics<<endl;
	cerr<<"xmin = "<<parametres.xmin<<" xmax = "<<parametres.xmax<<endl;
	cerr<<"ymin = "<<parametres.ymin<<" ymax = "<<parametres.ymax<<endl;
	//Initialisation a partir du fichier
	grid.writeGrid("grid.txt");
	return 0;
}


