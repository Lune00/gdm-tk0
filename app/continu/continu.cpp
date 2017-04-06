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
#include "vecteur.hpp"

struct Globals{
	int nx ;
	int ny ;
	int istart, iend, di;
	double xmin, xmax, ymin, ymax ;
	string metrics;
	string fichsim;
};


struct Champs{

	double * masse ;
	Vecteur * momentum ;
	Vecteur * vitesse ;
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

//Construction & initialisation de la grille a l'aide des globaux (lu)
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

//Initialisation uniquement ici des parametres globaux de l'analyse
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
				if(token=="Simu") is >> parametres.fichsim;
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
	cerr<<"nx = "<<parametres.nx<<" ny = "<<parametres.ny<<endl;
	cerr<<"xmin = "<<parametres.xmin<<" xmax = "<<parametres.xmax<<endl;
	cerr<<"ymin = "<<parametres.ymin<<" ymax = "<<parametres.ymax<<endl;
	Grid grid(parametres);
	grid.writeGrid("grid.txt");

	Simulation * mySimu = new Simulation();
	mySimu->read_data(parametres.fichsim.c_str());

	char nomFichier[100];

	for (unsigned int i = parametres.istart ; i != parametres.iend; i++){
		sprintf(nomFichier,"spl_nwk/spl_nwk_%.4d.his",i);
		cerr<<"****** Chargement : "<<nomFichier<<endl;
		mySimu->load_history(nomFichier);
		mySimu->algo()->algoFill();
	}

	return 0;
}


