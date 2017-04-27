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
	true_resolution_ = resolution_ * 0.5 ;

	for(unsigned int i = 0 ; i < nx_ ; i++)
	{
		for(unsigned int j=0 ; j < ny_;j++)
		{
			champ[ i * ny_ + j ] = 0. ;
		}
	}

}


Champ_Vectoriel::Champ_Vectoriel(unsigned int nx, unsigned int ny, string name, typechamps t,double resolution)
{
	nx_ = nx ;
	ny_ = ny ;
	name_ = name ;
	champx = new double [ nx_ * ny_ ];
	champy = new double [ nx_ * ny_ ];
	type_ = t ;
	resolution_ = resolution;
	//A mofifier ensuite: variance de l'exponentielle, mettre un parametres explicite
	true_resolution_ = resolution_ * 0.5 ;

	for(unsigned int i = 0 ; i < nx_ ; i++)
	{
		for(unsigned int j=0 ; j < ny_;j++)
		{
			champx[ i * ny_ + j ] = 0. ;
			champy[ i * ny_ + j ] = 0. ;
		}
	}

}
Champ_Scalaire::~Champ_Scalaire()
{
	cerr<<"Champ scalaire "<<this->name_<<" libéré."<<endl;
	delete [] champ ;
}

Champ_Vectoriel::~Champ_Vectoriel()
{
	cerr<<"Champ vectoriel "<<this->name_<<" libéré."<<endl;
	delete [] champx ;
	delete [] champy ;
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



void Champ_Vectoriel::calculMomentum(const Grid& grid,Sample& spl)
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
				champx[ i * ny_ + j ] += spl.body(id)->vx() * poids;
				champy[ i * ny_ + j ] += spl.body(id)->vy() * poids;
				poidstotal += poids;
			}

			if(poidstotal > epsilon_)
			{
				champx[ i * ny_ + j ] /= poidstotal ; 
				champy[ i * ny_ + j ] /= poidstotal ; 
			}
			else
			{
				champx[ i * ny_ + j ] = 0. ; 
				champy[ i * ny_ + j ] = 0. ; 
			}
		}
	}
}


void Champ_Vectoriel::calculVitesse(const Champ* masse,const Champ* momentum)
{
	for(unsigned int i = 0 ; i < nx_; i++)
	{
		for(unsigned int j=0 ; j < ny_ ; j++) 
		{
			champx[i*ny_+j] = momentum->getchamp_ij_x(i,j)/masse->getchamp_ij(i,j);
			champy[i*ny_+j] = momentum->getchamp_ij_x(i,j)/masse->getchamp_ij(i,j);
		}
	}
}


void Champ_Scalaire::writechamp(const Grid& grid)
{
	string filename = name_ + ".txt" ;
	ofstream champout (filename);
	for(int j = 0 ; j!= ny_  ; j++)
		for(int i = 0 ; i != nx_ ; i++)
		{
			{
				champout<<i+1<<" "<<j+1<<" "<<grid.getX(i,j)<<" "<<grid.getY(i,j)<<" "<<champ[ i * ny_ + j]<<endl;
			}
		}
	champout.close();

}

void Champ_Vectoriel::writechamp(const Grid& grid)
{
	string filename = name_ + ".txt" ;
	ofstream champout (filename);
	for(int j = 0 ; j!= ny_  ; j++)
		for(int i = 0 ; i != nx_ ; i++)
		{
			{
				champout<<i+1<<" "<<j+1<<" "<<grid.getX(i,j)<<" "<<grid.getY(i,j)<<" "<<champx[ i * ny_ + j]<<" "<<champy[ i * ny_ + j]<<endl;
			}
		}
	champout.close();

}
void Champ_Scalaire::drawchamp()
{
	//	string name = name_ + ".bmp" ;
	//	bitmap_image image( nx_ , ny_ );
	//
	//	if(!image)
	//	{
	//		cerr<<"Error - Failed to open: "<<name<<endl;
	//	}
	//
	//	double min = champ [0];
	//	double max = champ [0];
	//	for(unsigned int i = 0 ; i< nx_ ; i++)
	//	{
	//		for(unsigned int j = 0 ; j<ny_ ; j++)
	//		{
	//			min = min > champ [ i * ny_ + j ] ? champ [ i * ny_ + j ] : min; 
	//			max = max < champ [ i * ny_ + j ] ? champ [ i * ny_ + j ] : max; 
	//		}
	//	}
	//	cerr<<"min="<<min<<" max="<<max<<endl;
	//
	//	//On teste palette jet : entre 0 et 999
	//	for(unsigned int i = 0 ; i!=nx_; i++)
	//	{
	//		for(unsigned int j = 0 ; j!=ny_ ; j++)
	//		{
	//			double value =  champ[i * ny_ + j] ;
	//			int index = (value - min) * 999 / (max - min) ;
	//			image.set_pixel( i , j , jet_colormap[index] );
	//		}
	//	}
	//	image.save_image(name);
	//
}



//if( i == 100 && j == 100)
//{
//	ofstream test("pointp.txt");
//	for(map<int,double>::iterator it = local.begin(); it != local.end(); it++)
//	{
//		int id = it->first ;
//		test<<grid.readPoint(i,j).getX()<<" "<<grid.readPoint(i,j).getY()<<" "<<spl.body(id)->x()<<" "<<spl.body(id)->y()<<" "<<spl.body(id)->sizeVerlet()<<endl;
//	}
//	test.close();
//}
//
