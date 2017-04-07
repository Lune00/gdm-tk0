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

//Parametres globaux de l'analyse
//Differents types de champs calculables:

enum typechamps{t_masse,t_momentum};

class Config{

	private:

		int nx ;
		int ny ;
		int istart, iend, di;
		double xmin, xmax, ymin, ymax ;
		string metrics;
		string fichsim;
		bool calc_masse;

	public:
		Config();
		~Config(){};
		void init(ifstream&);
		void readMetrics(string);

		int getnx() {return nx ;}
		int getny() {return ny ;}
		int getistart() {return istart;}
		int getiend() {return iend;}
		int getdi() {return di;}

		double getxmin() {return xmin ;}
		double getymin() {return ymin ;}
		double getxmax() {return xmax ;}
		double getymax() {return ymax ;}
		string getfichsim() {return fichsim;}

		bool getcalcmasse() {return calc_masse;}
};

Config::Config()
{
	nx = 0 ;
	ny = 0 ;
	istart = 0 ;
	iend = 0 ;
	di = 0 ;
	xmin = 0. ;
	xmax = 0. ;
	ymin = 0. ;
	ymax = 0. ;
	metrics = "";
	fichsim = "Simu.sim";
	calc_masse = false ;
}

//Initialisation uniquement ici des parametres globaux de l'analyse
void Config::init(ifstream& is)
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
				if(token=="nx") is >>nx ; 
				if(token=="ny") is >>ny ; 
				if(token=="metrics") 
				{
					is >> token;
					readMetrics(token);
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
				if(token=="ini") is >>istart ; 
				if(token=="end") is >>iend ; 
				if(token=="di") is >>di; 
				if(token=="Simu") is >> fichsim;
				else if (token=="}") break;
				is >> token;
			}
		}
		if(token=="Champs{")
		{
			is >> token;
			while(is)
			{
				if (token=="masse") calc_masse = true ;
				if (token=="}") break;
				is >> token;
			}
		}

		is >> token;
	}
}

void Config::readMetrics(string file)
{
	ifstream is(file);
	string token;
	is >> token;
	while(is)
	{
		if(token=="xmin") is >> xmin;
		if(token=="xmax") is >> xmax;
		if(token=="ymin") is >> ymin;
		if(token=="ymax") is >> ymax;
		is >> token;
	}
}

class Champ
{

	protected:
		unsigned int nx_ ;
		unsigned int ny_ ;
		string name_ ;
		typechamps type_;//identifie le champ pour le calcul
	public:
		Champ();
		virtual ~Champ();
		virtual void calculMasse() {};
		string getname() {return name_ ; }
		typechamps gettype() {return type_ ; }
};

Champ::Champ(){}


Champ::~Champ(){}



class Champ_Scalaire : public Champ
{
	private:
	double * champ ;
	public :
	Champ_Scalaire(unsigned int,unsigned int, string,typechamps);
	~Champ_Scalaire();
	void calculMasse();

};

Champ_Scalaire::Champ_Scalaire(unsigned int nx, unsigned int ny, string name, typechamps t)
{
	nx_ = nx ;
	ny_ = ny ;
	name_ = name ;
	champ = new double [ nx_ * ny_ ];
	type_ = t ;
}

Champ_Scalaire::~Champ_Scalaire()
{
	cerr<<"Champ scalaire "<<this->name_<<" libéré."<<endl;
	delete [] champ ;
}

void Champ_Scalaire::calculMasse()
{

}
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

//Grille d'interpolation
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

//Construction & initialisation de la grille a l'aide des globaux (lu)
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

int main (int argc,char **argv)
{
	Config parametres;
	ifstream is(argv[1]);
	parametres.init(is) ;

	Grid grid(parametres);
	grid.writeGrid("grid.txt");

	Simulation * mySimu = new Simulation();
	mySimu->read_data(parametres.getfichsim().c_str());

	char nomFichier[100];

	for (unsigned int i = parametres.getistart() ; i != parametres.getiend(); i+=parametres.getdi()){
		sprintf(nomFichier,"spl_nwk/spl_nwk_%.4d.his",i);
		cerr<<"****** Chargement : "<<nomFichier<<endl;
		mySimu->load_history(nomFichier);
		mySimu->algo()->algoFill();

		//Calcul des champs
	}

	return 0;
}


