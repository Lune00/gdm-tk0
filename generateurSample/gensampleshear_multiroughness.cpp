#include<iomanip>
#include<cmath>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

class Particule{

  private:
    double x_, y_,r_;
    double rot_, vx_, vy_, vrot_;
    string type_;
    unsigned int group_;
  public:
    Particule();
    ~Particule();
    double getr() const {return r_;};
    void setr(double r) { r_ = r ;};
    double getx() const {return x_;};
    void setx(double x) { x_ = x ;};
    double gety() const {return  y_;};
    void sety(double y) { y_ = y ;};
    double getrot() const {return rot_;};
    double getvrot() const {return vrot_;};
    double getvx()  const {return vx_;};
    double getvy() const {return vy_;};
    std::string gettype() const {return type_;};
    double getgroup()  const {return group_;};
};


Particule::Particule()
{
  type_ = "disk"; 
  group_ = 0 ; 
  x_ = 0. ;
  y_ = 0. ;
  r_ = 0. ;
  vx_ = 0. ;
  vy_ = 0. ;
  rot_ = 0. ;
  vrot_ = 0. ;
};

Particule::~Particule() {};



double rand(double min,double max){
  return ((max-min)*((double)rand()/(double)RAND_MAX)+min);
}

double powerlaw(double r1, double r2){
  double y = rand(0.,1.);
  double tmp = (1./r2-1./r1)*y + 1./r1;
  return 1./tmp;
}


int main(){


  //Paramètres echantillon et de la rugosite de la paroi

  //Echantillon:

  unsigned int nfree = 25000 ;
  unsigned int ncouches = 150 ;
  double const r1 = 0.01 ;
  double const r2 = 2 * r1 ;
  double const rmean = 0.5 * ( r1 + r2 ) ;

  //Propriétés des plaques:

  //Plaque superieure:
  double const deltawS = 40;
  double const K_sup = 1.7 ;
  double const b = sin(deltawS*(M_PI/180.))/K_sup;
  double const R_sup = b/(1-b) ; // R=rparoi/rmean
  double const rparoi_sup = R_sup * rmean ;
  double const lw_sup =  K_sup * (2. * rparoi_sup) ;

  //Plaque inférieure:
  //Angle desire, en degre, K fixe, R ajuste
  double const deltaw = 22;
  double const K_inf = 1.7 ; // R=rparoi/rmean
  double const a = sin(deltaw*(M_PI/180.))/K_inf;
  double const R_inf = a/(1-a) ; // R=rparoi/rmean
  double const rparoi_inf = R_inf * rmean ;
  double const lw_inf =  K_inf * (2. * rparoi_inf) ;


  /* * * * * * * * * *  */

  unsigned int nslice = nfree / ncouches ;


  if(nfree == 0) {cout<<" Entrez un nombre de particules différent de 0."<<endl; return 0 ; }
  if(nslice ==0) {cout<<" Le paramètre nslice a une valeur non acceptable."<<endl; return 0; }
//  if(u < 0. || u * 2 * rparoi > 2 * r1) {cout<<"L'espacement entre particules de la paroi est trop grand (trous) ou trop faible (overlaping)."<<endl; return 0;}

  cout<<"--- Paroi superieure ---"<<endl;
  cout<<"Rugosite apparente R paroi : "<<R_sup<<endl;
  cout<<"Espacement entre particules de la paroi lw : "<< lw_sup <<endl;
  cout<<"Angle de rugosite paroi sup : "<<asin( K_sup * R_sup / ( 1. + R_sup) ) * 180. / M_PI <<" °"<<endl;

  cout<<"--- Paroi inferieure ---"<<endl;
  cout<<"Rugosite apparente R paroi : "<<R_inf<<endl;
  cout<<"Espacement entre particules de la paroi lw : "<< lw_inf <<endl;
  cout<<"Angle de rugosite paroi inf : "<<asin( K_inf * R_inf / ( 1. + R_inf) ) * 180. / M_PI <<" °"<<endl;

  std::vector<Particule> sample(nfree);
  std::vector<Particule> paroi1;
  std::vector<Particule> paroi2;

  cout<<"--- Echantillon ---"<<endl;
  cout<<"Nombre de particules libres : "<<sample.size()<<endl;
  cout<<"Rayon moyen des particules : "<<rmean<<endl;

  double xmin,xmax,ymin,ymax;

  xmin = 0. ;
  ymin = 0. ;
  // Initiatlisation des rayons:
  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  it->setr( powerlaw(r1,r2));
  }
  // Initialisation des positions sur la grille:
  unsigned int k = 0 ;
  unsigned int j = 0 ;
  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  double x, y;
	  double eps = rand(0.,1.) * r1 * 0.5   ;

	  if(it == sample.begin() || j == 0) 
	  { 
		  x = xmin ;
		  y = ymin + k * 2 * r2 + eps ;
	  }
	  else
	  {
		  x = (it - 1)->getx() + (it - 1)->getr() + it->getr() + eps ;
		  y = ymin + k * 2 * r2 + 0.5 * r1  ;
	  }

	  j++;

	  if( j == nslice ) { k++; j = 0 ; }

	  it->setx(x);
	  it->sety(y);
  }

  // Initialisation des parois et de leur position

  xmax = sample[0].getx();
  ymax = sample[0].gety();

  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  xmax = it->getx() > xmax ? xmax = it->getx() : xmax ;
	  ymax = it->gety() > ymax ? ymax = it->gety() : ymax ;
  }

  //Avoid overlaping 
  ymax += 3 * rmean ;
  ymin -= 3 * rmean ;
  xmin -= 3 * rmean ;
  xmax += 3 * rmean ;

  cout<<"xmax = "<<xmax<<endl;
  cout<<"ymax = "<<ymax<<endl;

  double lx = xmax - xmin ;

  unsigned int nparoi_sup = ceil( lx / ( lw_sup )) + 2;
  unsigned int nparoi_inf = ceil( lx / ( lw_inf )) + 5;

  cout<<"Nombre de particules composant la paroi sup : "<<nparoi_sup<<endl;
  cout<<"Nombre de particules composant la paroi inf : "<<nparoi_inf<<endl;

  paroi1.resize(nparoi_inf);
  paroi2.resize(nparoi_sup);

  // Initiatlisation des rayons des parois inf et sup:
  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
	  it->setr( rparoi_inf );
  }
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
	  it->setr( rparoi_sup );
  }

  // Paroi inferieure 
  k = 0 ;
  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
	  double x,y;
	  if( it == paroi1.begin() )
	  {
		  x = xmin + k * 2 * rparoi_inf ;
	  }
	  else
	  {
		  x = xmin + k * lw_inf ; 
	  }
	  y = ymin;
	  it->setx(x);
	  it->sety(y);
	  k++;
  }


  // Paori supérieure
  k = 0 ;
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
	  double x,y;
	  if( it == paroi2.begin() )
	  {
		  x = xmin + k * 2 * rparoi_sup ;
	  }
	  else
	  {
		  x = xmin + k * lw_sup ;
	  }

	  y = ymax;
	  it->setx(x);
	  it->sety(y);
	  k++;
  }
  //Ecrture du fichier packing0.spl

  ofstream myFile ("packing0.spl",ios::out);
  myFile<<"Sample{"<<endl;
  myFile<<"Cluster{"<<endl;
  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
	  myFile<<it->gettype()<<" "<<it->getgroup()<<" "<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }
  myFile<<"}"<<endl;
  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  myFile<<it->gettype()<<" "<<it->getgroup()<<" "<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }
  myFile<<"Cluster{"<<endl;
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
	  myFile<<it->gettype()<<" "<<it->getgroup()<<" "<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }
  myFile<<"}"<<endl;
  myFile<<"}"<<endl;
  myFile.close();

  ofstream plotsample("sample-gnuplot.txt");

  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
	  plotsample<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }

  for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
	  plotsample<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
	  plotsample<<it->getr()<<" "<<it->getx()<<" "<<it->gety()<<" "<<it->getrot()<<" "<<it->getvx()<<" "<<it->getvy()<<" "<<it->getvrot()<<endl;
  }
  plotsample.close();
  return 0;
}
