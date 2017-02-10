#ifndef _brazilian_A_h
#define _brazilian_A_h


#include "system_A.hpp"


#include <fstream>
//#include "dir.h"

#include "brazilian.hpp"

#include "body2d.hpp"
#include "inter2d.hpp"
//#include "heightProbe.hpp"
#include "circularProbe.hpp"
#include "rectangularProbe.hpp"
//#include "io.hpp"
//#include "system.hpp"
#include "solidfraction.hpp"
#include "anisotropy.hpp"
#include "stress.hpp"
#include "tensor.hpp"
#include "speedProfile.hpp"
#include "dataSet.hpp"
//#include "pointSet.hpp"
#include "NRsource.hpp"
#include "cluster_A.hpp"
#include "polyg.hpp"
#include "vertex.hpp"
#include "groupRelationData.hpp"


using namespace std;


class brazilian_A : public System_A
{

protected:

	//shearP_CD *sys;
	

	double time;

 	circularProbe totalProbe_;
	rectangularProbe prb_, probe, prb;
   
    rline * bottom_;
	rline * top_;
	rline * left_;
	rline * right_;
	GroupRelationData * grpRel_;
	unsigned int topId,bottomId,leftId,rightId;
	

//Donnees d'analyse du systeme
	unsigned int Nanalyze_;
	double sf_;
	double a_,oa_,an_,at_,al_;
	double da_,dn_,dt_,ds_,dl_;
	double qop_;
	double Z_,ratlers_;
	double fracdim_;
	double gaprmin,gaprmax,gapmoy_;
	double s1,s2;
	unsigned int Ngapsup,Ngap;
	double Vmoy_;
	double initY,initX,X0,Y0;
	double incR;
	
	double fnmoy_;
	double pressure_;

	double l0,h0;
	double defx,defy;
	double epsp,epsq;
	
	double zoom_;
	double ratio_;
	unsigned int divisionl;
	unsigned int divisionh;
	unsigned int nfracture;
	
	gdm::Tensor2x2  strain;
	
	vector< body2d *> rattlersB;
	vector< bool> rattlersVisu;
	
	//Pour l'affichage
	vector< bool > PGTM_;//greater than mean pressure in the system : allocated in granulostress2
	vector< bool > FSGTM_;//Force Sum greater than mean sum force in the system : allocated in forceMaxcorrelation
	//vector< bool > 

	body2d   * ngap1_, *ngap2_, *body_;
	inter2d	 * inter;
//ofstream data("Analyze.dat",ios::iout);
	bool initTime_;
	
	bool calcsf;
	bool calcFabric ;
	bool calcforcesA ;
	bool calcFabricPolyg;
	bool calcz ;
	bool calczp;
	bool calcfn;
	bool calcft;
//	bool calcSprofile;
//	bool calcSFprofile;
	bool calcglobalstress;
	bool calcgranulostress;
	bool normpdf;
	bool calcgap;
	bool calcdef;
	bool calcPtheta;
	bool calcFC;
	bool calczg;
	bool calcfracdim;
	bool granulo;
//	bool calcrfd;
	bool calcgranulopdf;
	
	bool removeR;
	bool growR;
	bool removeVisu;
	
	bool pressure_defined;
	bool sumforce_defined;
	
	// modif du 17/03/09 
	bool calClust;
	//bool calcPtheta_Clust;
	cluster_A* cla_;
	
	bool displaySample;
	bool displayBranch;
	bool displayForce;
	bool calcompact;
	bool ratio;
	bool calfracture;
	
	unsigned int NFabric;
	int NbinFN;
	int NbinFA;
	int NbinFT;
	int NbinPT;
	int Nprb_;
	int Nbingranulo;
	int mobperiod,mobwidth;
	int perF,wF;
	int perG,wG;
	int Nd,Np;
	int Nq_;
	long int numFile_;
	unsigned int NbinClus,mobperiodClus,mobwidthClus;


public:
	
    void initAnalyse( );
	void analyse(double, unsigned int, unsigned int);
	void read_parameters(istream &);
	void plugRef();
	
	//void plug( shearP_CD * sys_in) {sys=sys_in;}

	~brazilian_A() { }
	brazilian_A(){calcsf=calcFabric=calcz=calczp=calcfn=calcft=calcFC=calczg=
		calcglobalstress=calcgranulostress=normpdf=calcforcesA=calcPtheta=calcfracdim=granulo=
		removeR=growR=calcgranulopdf=pressure_defined=sumforce_defined=removeVisu=calcFabricPolyg=displaySample=displayBranch=displayForce=calfracture=false;
	calcgap=calcdef=initTime_=true;
	fnmoy_=1.;
	zoom_ = 1.;
	numFile_ = 0;
	nfracture = 0;
	}
	
	brazilian_A(brazilian *sys_a) { sys_=sys_a; calcsf=calcFabric=calcz=calczp=calcfn=calcft=calcFC=calczg=
		calcglobalstress=calcgranulostress=normpdf=calcforcesA=calcPtheta=calcfracdim=granulo=
		removeR=growR=calcgranulopdf=pressure_defined=sumforce_defined=removeVisu=calcFabricPolyg=displaySample=displayBranch=displayForce=calfracture=false;
	calcgap=calcdef=initTime_=true;
	fnmoy_=1.;
	zoom_ = 1.;
	numFile_ = 0;
	}

	
	void allFalse() {calcsf=calcFabric=calcz=calczp=calcfn=calcft=calcFC=calczg=
		calcglobalstress=calcgranulostress=normpdf=calcforcesA=calcPtheta=calcfracdim=granulo=
		removeR=growR=calcgranulopdf=pressure_defined=sumforce_defined=removeVisu=calcFabricPolyg=displaySample=displayBranch=displayForce=calfracture=false;
	calcgap=calcdef=initTime_=true;
	fnmoy_=1.;
	zoom_ = 1.;
	numFile_ = 0;
	}

	circularProbe & totalProbe()  {return  totalProbe_;}
	rectangularProbe & rProbe()  {return  prb_;}

	unsigned int & Nanalyze()          { return Nanalyze_; }
	unsigned int   Nanalyze()    const { return Nanalyze_; }

	double & sf()          { return sf_; }
	double   sf()    const { return sf_; }

	double & a()          { return a_; }
	double   a()    const { return a_; }

	double & an()          { return an_; }
	double   an()    const { return an_; }

	double & at()          { return at_; }
	double   at()    const { return at_; }
	
	double & al()          { return al_; }
	double   al()    const { return al_; }

	double & z()          { return Z_; }
	double   z()    const { return Z_; }
	
	bool & initTime()   {return initTime_;}
	bool   initTime()   const {return initTime_;}
	

//	body2d *  Partref()  const {return partref;}
//	body2d * &  Partref()  {return partref;}
	body2d *  ngap1()  const {return ngap1_;}
	body2d *  ngap2()  const {return ngap2_;}
	
//le 06/03/09
//accesseur de la classe cluster_A
	cluster_A * clusterA() const {return cla_; }

//Analysis functions
	void Gap();
	void def();
	void SF();
	void A();
	void forces_A( int );
	unsigned int Z(bool,unsigned int,bool);
	void globalStress();
	
	//void profiles(bool , bool);
	int  granulopdf(unsigned int,bool,int,bool,unsigned int, unsigned int );
	int  pdfforce(bool,int,bool,unsigned int, unsigned int);//fn?  nbin normalized?
	void Ptheta( unsigned int,unsigned int,unsigned int);
	void forcesMaxCorrelation();
	void granuloSpeed(unsigned int );
	
	void polygAnisotropy(unsigned int);
	
	/*void fractalDimension( );
	void fractalDensity();
	void BboxCounting();
	void NboxCounting();
	void BmultifractalBC( );
	void NmultifractalBC( );
	*/
	
	//void rfd( unsigned int, unsigned int);
	
	void removeBody( body2d* );
	void removeRattlers( );
	void reduceRattlers();
	void growRattlers(double );

	unsigned int granuloStress( unsigned int );//Nbin
	unsigned int granuloStress2( unsigned int ,unsigned int,unsigned int );//Nevent per size classes
	unsigned int granuloStress3( unsigned int );//Nbin sliding histo
	
	void writePS(const char * fname);
	void writePS2(const char * fname);
	void display_data(char * fname);
	
	void cluster();
	bool Inter (polyg *p, polyg *q);
	polyg * Intersection (polyg *p, polyg *q);
	//double compactness(double &x);
	double compactness(double &h, double &l);
	double compactness(double &x, double &y, double &h, double &l);
	double compactness_disk(double &x, double &y, double &h, double &l);
	void fracture();
};

#endif 



