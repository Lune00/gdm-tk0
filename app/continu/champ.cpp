#include"champ.hpp"

using namespace std;

Champ::Champ(){ epsilon_ = 0.0000001; }

Champ::~Champ(){}

//Renvoie la ponderation exponentielle a partir de la distance au point et de la resolution (true_resolution)
double Champ::ponderation(double d)
{
	return exp( (-d/true_resolution_) * (-d/true_resolution_) );
}

Champ_Scalaire::Champ_Scalaire(unsigned int nx, unsigned int ny, string name, typechamps t,double resolution)
{
	nx_ = nx ;
	ny_ = ny ;
	name_ = name ;
	champ = new double [ nx_ * ny_ ];
	type_ = t ;
	resolution_ = resolution;
	//A mofifier ensuite: variance de l'exponentielle, mettre un parametres explicite
	true_resolution_ = resolution_ * 0.8;

	for(unsigned int i = 0 ; i < nx_ ; i++)
	{
		for(unsigned int j=0 ; j < ny_;j++)
		{
			champ[ i * ny_ + j ] = 0. ;
		}
	}

}


Champ_Scalaire::~Champ_Scalaire()
{
	cerr<<"Champ scalaire "<<this->name_<<" libéré."<<endl;
	delete [] champ ;
}

void Champ_Scalaire::calculMasse(const Grid& grid,Sample& spl)
{
	for(unsigned int i = 0 ; i < nx_; i++)
	{
		for(unsigned int j=0 ; j < ny_ ; j++) 
		{
			double poidstotal = 0. ;
			map<int,double> local = grid.readPoint(i,j).getparticules();

			for(map<int,double>::iterator it = local.begin(); it != local.end(); it++)
			{
				int id = it->first ;
				double d = it->second;
				double poids = ponderation(d);
				champ[ i * ny_ + j ] += spl.body(id)->mass() * poids;
				poidstotal += poids;
			}

			if(poidstotal > epsilon_) champ[ i * ny_ + j ] /= poidstotal ; 
			else
				champ[ i * ny_ + j ] = 0. ; 
		}
	}
}

void Champ_Scalaire::writechamp(const Grid& grid)
{
	string filename = name_ + ".txt" ;
	ofstream champout (filename);
	for(int j = 0 ; j != ny_ ; j++)
	{
		for(int i = 0 ; i!= nx_  ; i++)
		{
			champout<<grid.getX(i,j)<<" "<<grid.getY(i,j)<<" "<<champ[ i * ny_ + j]<<endl;
		}
	}

}
