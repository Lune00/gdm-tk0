#include<iomanip>
#include<random>
#include<limits>
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
    void setgroup(unsigned int id) {group_ = id;};
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

  //user
  //Paramètres echantillon et de la rugosite de la paroi
  double const r1 = 0.01 ;
  double const r2 = 3 * r1 ;
  double const rmean = 0.5 * ( r1 + r2 ) ;
  double const R = 0.3 ; // R=rparoi/rmean
  double const u = 0.2 ; // lw = 2*rparoi + u * paroi
  double const rparoi = R * rmean ;
  double const lw = 2. * rparoi ;// (1. + u ) * ( 2 * r1 );

  unsigned int nfree = 500;
  unsigned int nslice = 50 ;

  //GroupData:
  double masse_vol = 2600.0 ;

  //GroupRelationData:
  double const ah = 0. ;
  double const en = 0. ;
  double const et = 0. ;
  double const rs = 0. ;
  // nombre de classes de particules differentes uniquement pour les particules libres  (assigner des valeurs de mu, cohesion differentes par ex)
  //si ngroup = 1 alors il y a 2 groupes : particules libres et parois
  unsigned int ngroup = 3 ;
  //Par defaut les partois ont un groupe different
  unsigned int groupparoi = ngroup ; // par defaut

  //Assigne a chaque paire group1-group2 une valeur de mu_s (frottement sec) suivant une distribution normale N(mu_s_moyen,mu_s_variance)
  double mu_s_moyen = 0.3 ;
  double mu_s_variance = 0.1 ;
  // le frottement par defaut (true) entre grains libres et parois est identique pour tous les groupes, sinon il suit la loi des combinaisons possibles
  // si mu_s_paroi_any est vrai alors le frottement entre (ngroup+1 <-> n) est egal a mu_s_moyen  
  bool mu_s_paroi_any = true ;
  //Ignore tous les groupes:
  bool uniform_mu_s_for_all = false ;

  int seed = 3 ;

  //Random generator:
  srand(seed);
  std::default_random_engine generator;

  std::normal_distribution<double> distribution(mu_s_moyen,mu_s_variance);
  if(ngroup == 0 ) {cout<<"ngroup ne peut pas être egal à 0 ! Il doit y avoir au moins un groupe."<<endl;return 0;}
  if(nfree == 0) {cout<<" Entrez un nombre de particules différent de 0."<<endl; return 0 ; }
  if(nslice ==0) {cout<<" Le paramètre nslice a une valeur non acceptable."<<endl; return 0; }
  if(u < 0. || u * 2 * rparoi > 2 * r1) {cout<<"L'espacement entre particules de la paroi est trop grand (trous) ou trop faible (overlaping)."<<endl; return 0;}
  if(ngroup > 50 ) {cout<< "Attention, le nombre de groupe est élevé, il y a "<<ngroup*ngroup<<" possibilités d'interaction ! "<<endl;}
  cout<<"Rugosite apparente R : "<<R<<endl;
  cout<<"Espacement entre particules de la paroi lw : "<< lw <<endl;
  cout<<"Angle de rugosite : "<<asin( lw / ( (2 * rmean) * ( 1 + R ) ) ) * 180. / M_PI <<" °"<<endl;

  std::vector<Particule> sample(nfree);
  std::vector<Particule> paroi1;
  std::vector<Particule> paroi2;
  cout<<"Nombre de particules libres : "<<sample.size()<<endl;
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
  ymax += 2 * rmean ;
  ymin -= 2 * rmean ;
  xmin -= 2 * rmean ;
  xmax += 2 * rmean ;

  cout<<"xmax = "<<xmax<<endl;
  cout<<"ymax = "<<ymax<<endl;

  double lx = xmax - xmin ;
  unsigned int nparoi = ceil( lx / ( lw )) + 2;

  cout<<"Nombre de particules composant la paroi : "<<nparoi<<endl;

  paroi1.resize(nparoi);
  paroi2.resize(nparoi);

  // Initiatlisation des rayons des parois sup et inf:
  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
    it->setr( rparoi );
  }
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
    it->setr( rparoi );
  }

  // Paroi inferieure 
  k = 0 ;
  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
    double x,y;
    if( it == paroi1.begin() )
    {
      x = xmin + k * 2 * rparoi ;
    }
    else
    {
      x = xmin + k * lw ; // ( 2 * rparoi + lw ) ;
    }
    y = ymin;
    it->setx(x);
    it->sety(y);
    k++;
  }
  k = 0 ;
  // Paori supérieure
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
    double x,y;
    if( it == paroi2.begin() )
    {
      x = xmin + k * 2 * rparoi ;
    }
    else
    {
      x = xmin + k * lw ; //( 2 * rparoi + u ) ;
    }

    y = ymax;
    it->setx(x);
    it->sety(y);
    k++;
  }
  //Assignation des groupes de manière aléatoire pour les particules libres:

  if( ngroup != 0 ) { 
    for(std::vector<Particule>::iterator it = sample.begin() ; it!= sample.end(); it++){
      unsigned int group = rand () % ngroup + 1;
      it->setgroup(group);
    }

  }

  //Assignation du groupe aux particules des  parois
  for(std::vector<Particule>::iterator it = paroi1.begin() ; it!= paroi1.end(); it++){
    it->setgroup(groupparoi);
  }
  for(std::vector<Particule>::iterator it = paroi2.begin() ; it!= paroi2.end(); it++){
    it->setgroup(groupparoi);
  }
  //Generation d'un fichier group.ini pour les relations entre groupes a inclure dans Simu.sim
  ofstream groupIni("group.ini",ios::out);

  groupIni <<"GroupData{"<<endl;
  groupIni <<"ngrp "<<ngroup+1<<endl;
  groupIni <<"parameter density "<<endl;
  groupIni <<"setall density "<<masse_vol<<endl;
  groupIni <<"}"<<endl;

  groupIni<<"GroupRelationData{"<<endl;

  groupIni <<"ngrp "<<ngroup+1<<endl;

  groupIni <<"parameter ah "<<endl;
  groupIni <<"setall ah "<<ah<<endl;

  groupIni <<"parameter en "<<endl;
  groupIni <<"setall en "<<en<<endl;

  groupIni <<"parameter et "<<endl;
  groupIni <<"setall et "<<et<<endl;

  groupIni <<"parameter rs "<<endl;
  groupIni <<"setall rs "<<rs<<endl;

  groupIni <<"parameter mu"<<endl;
  if(uniform_mu_s_for_all) {
    groupIni<<"setall mu "<<mu_s_moyen<<endl;
  }
  else
  {
    for (unsigned int i = 0 ; i != ngroup + 1 ; i++){
      for (unsigned int j = i ; j != ngroup + 1 ; j++){

	double mu = distribution(generator);

	if( i == groupparoi || j == groupparoi  )
	{
	  if (mu_s_paroi_any){
	    groupIni<<"set mu "<<i<<" "<<j<<" "<<mu_s_moyen<<endl;
	  }
	  else{ 
	    groupIni<<"set mu "<<i<<" "<<j<<" "<<mu<<endl;
	  }
	}
	else
	{
	  groupIni<<"set mu "<<i<<" "<<j<<" "<<mu<<endl;
	}
      }
    }
  }




  groupIni<<"}"<<endl;
  groupIni.close();
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
