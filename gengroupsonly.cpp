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


int main() {

  int seed = 3 ;

  //Random generator:
  srand(seed);
  std::default_random_engine generator;

  // nombre de classes de particules differentes uniquement pour les particules libres  (assigner des valeurs de mu, cohesion differentes par ex)
  //si ngroup = 1 alors il y a 2 groupes : particules libres et parois
  unsigned int ngroup = 3 ;
  //Par defaut les partois ont un groupe different
  unsigned int groupparoi = ngroup ; // par defaut

  //GroupData : 
  double masse_vol = 2600.0 ;

  //GroupRelationData:
  double const ah = 0. ;
  double const en = 0. ;
  double const et = 0. ;
  double const rs = 0. ;

  //Assigne a chaque paire group1-group2 une valeur de mu_s (frottement sec) suivant une distribution normale N(mu_s_moyen,mu_s_variance)
  double mu_s_moyen = 0.3 ;
  double mu_s_variance = 0.1 ;
  // le frottement par defaut (true) entre grains libres et parois est identique pour tous les groupes, sinon il suit la loi des combinaisons possibles
  // si mu_s_paroi_any est vrai alors le frottement entre (ngroup+1 <-> n) est egal a mu_s_moyen  
  bool mu_s_paroi_any = true ;
  //Ignore tous les groupes:
  bool uniform_mu_s_for_all = false ;


  std::normal_distribution<double> distribution(mu_s_moyen,mu_s_variance);

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
}
