#include "grid.hpp"

int Grid::check(Config parametres)
{
	if(parametres.getl() < sqrt(dx_ * dx_ + dy_ * dy_ ) * 0.5 )
	{
		cerr<<"@initmotif : La resolution est trop faible par rapport au pas de la grille - recouvrement de l'espace partiel - augmenter spatial "<<endl;
		return 1;
	}
	else
		return 0;
}

Grid::Grid(double dx,double dy,double l, int mult)
{
	unsigned int nx = round(l / dx);
	unsigned int ny = round(l / dy);
	dx_ = dx ;
	dy_ = dy ;
	xmin_ = 0. ;
	ymin_ = 0. ;
	xmax_ = mult * nx ;
	ymax_ = mult * ny ;
	cerr<<"nx = "<<nx<<" - ny = "<<ny<<endl;
}

//renvoie vrai sur les deux points se recouvrent via la resolution spatiale
bool Grid::recouvrement(Point a, Point b) 
{
	double xa = a.getX();
	double ya = a.getY();
	double xb = b.getX();
	double yb = b.getY();

	double dist = sqrt( (xb - xa) * (xb - xa) + (yb - ya) * (yb - ya) ) ;
	if ( dist < resolution_ ) return true;
	else
		return false;
}


Point Grid::getPoint(int i,int j)
{
	if ( i >= 0 && i <= nx_ && j >=0 && j<= ny_) return array_[ i * ny_ + j ] ; 
	else
	{
		//Gerer exception:
		return array_[ i * ny_ + j ] ;
	}
}

void Grid::initmotif(Config& parametres)
{
	//On cree le motif a partir d'une grille tmp , que l'on appliquera a tous les domaines
	//Le motif nous renvoie les cordonnées relatives au point considéré (en bas a gauche d'un domaine auquel appartient la particule) de tous les points auquels cette particule apporte une contribution

	int size = 10 ;

	Grid *  tmp = new Grid(dx_,dy_,resolution_,size);

	tmp->setcoordinates();
	//On prend le point au centre de la grille:
	int iref = size / 2 ;
	int jref = size / 2 ;
	//On commence par stocker tous les points qui recouvrent ce point x
	//   o-------o
	//   |       |
	//   x-------o

	for (unsigned int j = 0 ; j!=tmp->ny_ ; j++)
	{
		for(unsigned int i = 0 ; i != tmp->nx_ ; i++)
		{
			cerr<<i<<" "<<j<<endl;
			if(tmp->recouvrement(tmp->getPoint(iref,jref),tmp->getPoint(i,j)))
			{
				cerr<<"Ca recouvre !"<<endl;
			}
		}
	}
	cerr<<"End initmotif"<<endl;
	delete tmp ;
}

void Grid::setcoordinates()
{
	//Set coordinates:
	double x = xmin_;
	double y = ymin_;

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

Grid::Grid(Config parametres)
{
	nx_ = parametres.getnx() ;
	ny_ = parametres.getny() ;

	xmin_ = parametres.getxmin() ;
	xmax_ = parametres.getxmax() ;
	ymin_ = parametres.getymin() ;
	ymax_ = parametres.getymax() ;
	resolution_ = parametres.getl();

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

	setcoordinates();
	initmotif(parametres);
	cerr<<"End Grid"<<endl;
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
