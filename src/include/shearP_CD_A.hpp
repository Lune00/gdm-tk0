#ifndef _shearP_CD_A_h
#define _shearP_CD_A_h


#include "system_A.hpp"


#include <fstream>
//#include "dir.h"


#include "shearP_CD.hpp"

#include "connectivity_A.hpp"


#include "body2d.hpp"
#include "inter2d.hpp"
#include "heightProbe.hpp"
#include "rectangularProbe.hpp"
#include "solidfraction.hpp"
#include "anisotropy.hpp"
#include "stress.hpp"
#include "tensor.hpp"
#include "speedProfile.hpp"
#include "dataSet.hpp"
//#include "pointSet.hpp"
#include "NRsource.hpp"

using namespace std;

class shearP_CD_A : public System_A
{

protected:

	//shearP_CD *sys;

	double time;

	//connectivity_A ca_;
	
//Donnees d'analyse du systeme
	unsigned int Nanalyze_;
	heightProbe totalProbe_;
	rectangularProbe totalProbe_R;
	double sf_;
	double a_,oa_,an_,at_,al_,anc_,atc_;
	double da_,dn_,dt_,ds_,dl_,dnc_,dtc_;
	double qop_,pressure_,q_;
	double Z_,ratlers_;
	double fracdim_;
	double gaprmin,gaprmax,gapmoy_;
	double s1,s2;
	unsigned int Ngapsup,Ngap;
	double Vmoy_;
	double initY,initX,X0,Y0;
	double incR;
	double X00,Y00;
	double Vinit_,dvov_;
	

	
	gdm::Tensor2x2  strain;
	double epsp,epsq;
	double epsxy;
	
	
	vector< body2d *> rattlersB;
	vector< bool> rattlersVisu;
	
	//Pour l'affichage
	vector< bool > PGTM_;//greater than mean pressure in the system : allocated in granulostress2
	vector< bool > FSGTM_;//Force Sum greater than mean sum force in the system : allocated in forceMaxcorrelation
	//vector< bool > 

	body2d * partref, * ngap1_, *ngap2_;
	unsigned int partrefId;
	//body2d * origindef;
	double w0;//h0;

	double fnmoy_;
	double lmoy_;
	
	bool calcsf;
	bool calcFabric ;
	bool plotFA;
	bool calcforcesA ;
	bool calcforcesAL ;
	bool calcforcesAC ;
	bool calcz ;
	bool calczp;
	
	bool calcfn;
	bool calcft;
	bool calcl;
	
	bool calcSprofile;
	bool calcSFprofile;
	bool calcglobalstress;
	bool calcpartialLengthstress;
	bool calcpartiaNormalForcestress;
	bool calcgranulostress;
	bool normpdf;
	bool calcgap;
	bool calcdef;
	bool calcPtheta;
	bool calcFC;
	bool calczg;
	bool calcfracdim;
	bool granulo;
	bool calcrfd;
	bool calcgranulopdf;
	bool exportDistribution;
	bool calcinout;
	//bool filter;
	
	bool removeR;
	bool growR;
	bool removeVisu;
	bool visu;
	
	bool pressure_defined;
	bool sumforce_defined;

	unsigned int Nbincor_;
	unsigned int Npointps;
	int NbinFN;
	int NbinL;
	int NbinFA;
	int NbinFT;
	int NbinPT;
	int Nprb_;
	int Nbingranulo;
	int mobperiod,mobwidth;
	int perF,wF;
	int perL,wL;
	int perG,wG;
	int Nd,Np;
	int Nq_;
	int Ninout;
	
	
	///////////
    bool ContactMesh_;
	bool displaySample;
	bool displayForce;
	double zoom_;

public:
	
    void initAnalyse( );
	void analyse(double, unsigned int, unsigned int);
	void read_parameters(istream &);
	void plugRef();

	
	//void plug( shearP_CD * sys_in) {sys=sys_in;}

	void allFalse() {calcsf=calcFabric=calcz=calczp=calcfn=calcft=calcSprofile=calcFC=calczg=
		calcSFprofile=calcglobalstress=calcgranulostress=normpdf=calcforcesA=calcPtheta=calcfracdim=granulo=calcrfd=
		removeR=growR=calcgranulopdf=pressure_defined=sumforce_defined=removeVisu=exportDistribution=calcl=visu=
		calcpartialLengthstress=calcpartiaNormalForcestress=displaySample=displayForce=ContactMesh_=false;
		calcgap=calcdef=true;
	fnmoy_=1.;
	Nanalyze_ = 0;
	zoom_ = 1;}
	
	~shearP_CD_A() { this->allFalse();	}
	shearP_CD_A(){this->allFalse();}
	
	shearP_CD_A(shearP_CD *sys_a) { this->allFalse();}

	


	heightProbe & totalProbe()  {return  totalProbe_;}

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
	
	

	body2d *  Partref()  const {return partref;}
	body2d * &  Partref()  {return partref;}
	body2d *  ngap1()  const {return ngap1_;}
	body2d *  ngap2()  const {return ngap2_;}

//Analysis functions
	void Gap();
	void def();
	void SF();
	void A();
	void forces_A( int,bool );
	void forces_AL( int,bool );
	void forces_AC( int,bool );
	
	unsigned int Z(bool,unsigned int,bool);
	void globalStress();
	void partialLengthStress(unsigned int);
	void partialNormalForceStress(unsigned int);

	void profiles(bool , bool);
	int  pdfforce(bool,int,bool,unsigned int, unsigned int);//fn?  nbin normalized?
	int  pdflength(int,bool,unsigned int, unsigned int);//  nbin normalized?
	void normalForceInOut(int );
	
	void Ptheta( unsigned int,unsigned int,unsigned int);
	void granuloSpeed(unsigned int );
	
	void removeBody( body2d* );
	void removeRattlers( );
	void reduceRattlers();
	void growRattlers(double );

	unsigned int granuloStress( unsigned int );//Nbin
	unsigned int granuloStress2( unsigned int ,unsigned int,unsigned int );//Nevent per size classes
	unsigned int granuloStress3( unsigned int );//Nbin sliding histo
	void followparticles();	
	void writePS(const char * fname);
	//void filterGap();



};

#endif // _shearP_CD_A_h



