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
	nx_ = mult * ceil(l / dx);
        ny_ = mult * ceil(l / dy);
	cerr<<nx_<<" - "<<ny_<<endl;
	dx_ = dx ;
	dy_ = dy ;
	xmin_ = 0. ;
	ymin_ = 0. ;
	xmax_ = mult * nx_ ;
	ymax_ = mult * ny_ ;
	resolution_ = l ;
	array_ = new Point [ nx_ * ny_ ];
	if(ny_ == 0 || nx_ == 0 ) cerr<<"@Grid::Grid(dx,dy,l,mult) grille temporaire - dimension nulle"<<endl;
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

bool sortv(const voisin& a,const voisin& b)
{
	return a.i_ < b.i_;
}



Point Grid::getPoint(int i,int j)
{
	if ( i >= 0 && i < nx_ && j >=0 && j< ny_) return array_[ i * ny_ + j ] ; 
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

	cerr<<"Initialisation du motif : taille "<<motif_.size()<<endl;
	//double dx = 1. ;
	//double dy = 2. ;
	//double resolution = 0.71;
	//if(resolution< sqrt(dx * dx + dy * dy ) * 0.5 ) cerr<<"RES TROP FAIBLE!"<<endl;
	Grid *  tmp = new Grid(dx_,dy_,resolution_,size);

	tmp->setcoordinates();
	tmp->writeGrid("tmp.txt");
	//On prend le point au centre de la grille:
	int iref = size / 2 ;
	int jref = size / 2 ;
	//On commence par stocker tous les points qui recouvrent ce point x
	//   o-------o
	//   |       |
	//   x-------o
	

	for(unsigned int l=jref; l!=jref + 2 ;l++)
	{
		for (unsigned int k = iref ; k != iref + 2 ; k++)
		{

			for (unsigned int j = 0 ; j!=tmp->ny_ ; j++)
			{
				for(unsigned int i = 0 ; i != tmp->nx_ ; i++)
				{
					if( i == k && j == l ) continue;

					if(tmp->recouvrement(tmp->getPoint(k,l),tmp->getPoint(i,j)))
					{
						voisin v;
						v.i_ = i - iref  ;
						v.j_ = j - jref  ;
						motif_.insert(v);
						if(!motif_.insert(v).second) {

							cerr<<"motif = { "<<endl;
							for(std::set<voisin>::iterator it = motif_.begin() ; it != motif_.end();it++)
							{	
								cerr<<"("<<(*it).i_<<","<<(*it).j_<<")"<<endl;
							}
							cerr<<"}"<<endl;

							cerr<<"Failed to insert ("<<v.i_<<","<<v.j_<<")"<<endl;
						}

					}
				}
			}
		}
	}

	//Par précaution on met toujours les 4 points du domaine de reference a tester:
	for(unsigned int l=jref; l!=jref + 2 ;l++)
	{
		for (unsigned int k = iref ; k != iref + 2 ; k++)
		{
			voisin v ;
			v.i_ = k - iref  ;
			v.j_ = l - jref  ;
			motif_.insert(v);
			if(!motif_.insert(v).second) {

				cerr<<"motif = { "<<endl;
				for(std::set<voisin>::iterator it = motif_.begin() ; it != motif_.end();it++)
				{	
					cerr<<"("<<(*it).i_<<","<<(*it).j_<<")"<<endl;
				}
				cerr<<"}"<<endl;

				cerr<<"Failed to insert ("<<v.i_<<","<<v.j_<<")"<<endl;
			}
		}
	}
	cerr<<"Taille du motif : "<<motif_.size()<<endl;
	ofstream imotif("motif_indices.txt");
	for(std::set<voisin>::iterator it = motif_.begin() ; it != motif_.end();it++)
	{	
		imotif<<(*it).i_<<" "<<(*it).j_<<" "<<(*it).i_*dx_+xmin_<<" "<<(*it).j_+ymin_<<endl;
	}
	imotif.close();

	delete tmp ;
}

//pas necessaire
bool Grid::motifinit()
{
	if(motif_.size()==0)
	{
		cerr<<"Erreur : motif ne contient aucun point d'intégration !"<<endl;
		return false;
	}
	else
		return true;
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
			array_[ i * ny_ + j ].seti(i);
			array_[ i * ny_ + j ].setj(j);
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


Point& Grid::returnPoint(int i,int j)
{
	return array_[ i * ny_ + j];
}

void Grid::writeGrid(string filename)
{
	ofstream gridout (filename,ios::out);
	for(int j = 0 ; j != ny_ ; j++)
	{
		for(int i = 0 ; i!= nx_ - 1 ; i++)
		{
			gridout<<getX(i,j)<<" "<<getY(i,j)<<" "<<dx_<<" 0."<<" "<<resolution_<<endl;
		}
	}

	for(int i = 0 ; i != nx_ ; i++)
	{
		for(int j = 0 ; j!= ny_ - 1 ; j++)
		{
			gridout<<getX(i,j)<<" "<<getY(i,j)<<" 0. "<<dy_<<" "<<resolution_<<endl; 
		}
	}

	gridout.close();
}	

//On stock les particules sur chaque point
void Grid::repartition(Sample& spl)
{
	unsigned int N = spl.lbody().size();
	cerr<<"Nombre de particules a stocker : "<< N <<endl;

	ofstream part("particules.txt");
	//On commence par parcourir les particules et a reperer leur point de ref
	for(unsigned int k = 1000 ; k != 1001 ; k++)
	{
		double x = spl.body(k)->x() ; 
		double y = spl.body(k)->y() ; 
		double r =  spl.body(k)->sizeVerlet() ; 
		unsigned int i = floor( (x-xmin_)/dx_);
		unsigned int j = floor( (y-ymin_)/dy_);
		part<<i * dx_  <<" "<<j* dy_ <<" "<<x -xmin_<<" "<<y- ymin_<<" "<<r<<endl;
		// (i,j) point de reference : on applique le motif a partir de ce point
		// Fonction qui prend la particule et la grille
		updatePoints(returnPoint(i,j),spl.body(k));
	}

	part.close();
}


double Grid::getdistance(Point p , body2d* b)
{
	return sqrt( (b->x() - p.getX())*(b->x() - p.getX()) + (b->y() - p.getY())*(b->y() - p.getY()) ) ;
}

bool Grid::belongtopoint(Point p , body2d* b)
{
	return (getdistance(p,b) < resolution_) ;
}

bool Grid::out(int i,int j)
{
	if( i < 0 || i > nx_ || j < 0 || j > ny_) return true;
	else
		return false;
}
void Grid::updatePoints(Point& ref,body2d* p)
{

	int iref = ref.geti() ;
	int jref = ref.getj() ;

	for(std::set<voisin>::iterator it = motif_.begin(); it != motif_.end() ; it++)
	{
		//Check si le motif sort de la grille:
		int i = iref + it->i_;
		int j = jref + it->j_;
		if(out(i,j)) continue;
		if( belongtopoint(returnPoint(i,j),p)) 
		{
			cerr<<"La particule "<<p->id()<<" appartient au point "<<i<<" "<<j<<endl;

		}


	}


}

