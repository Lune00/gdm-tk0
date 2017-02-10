#ifndef _cluster_A_hpp
#define _cluster_A_hpp

#include "system.hpp"
#include "interdof.hpp"
#include "probe.hpp"
#include "dataSet.hpp"
#include "rectangularProbe.hpp"
#include "tensor.hpp"

//Fonctions d'analyse statistique sur un groupe de cluster:
// topologie des contacts
// contrainte par cluster


//Question : quelles forces considérer : brute ou projetée sur le vect branche ??
//using namespace gdm;
class cluster_A

{
		vector<interdof*> linterdof_;
		System * sys_;
		unsigned int Ndof_;
		vector<bool> activatedDof;
		interdof * interdofk;
		dof * i_;
		dof * j_;
		double sigma_,ac_,an_,at_;
		double acs_,acns_,sigmaSimple_,sigmaNsimple_;
		vector<double>MeanRadial_;
		vector<double>MeanOrtho_;
		double Z;
		double Xratlers_;
		//double Xforce_simple_,Xforce_double_,Xforce_triple_;

				
public:
	
	cluster_A(System * sys/*, rectangularProbe & prb*/):sys_(sys)/*,prb_(prb)*/{sigma_ = ac_ = an_ = at_ = Z = acs_ = acns_ = Xratlers_ = 0;}
	void buildListinterdof( Probe &);
	vector<interdof*>& linterdof() {return linterdof_;}
	
	double sigma() {return sigma_;}
	double sigmaSimple()	{return sigmaSimple_;}
	double sigmaNsimple()	{return sigmaNsimple_;}
	
	double ac() {return ac_;}
	double an() {return an_;}
	double at() {return at_;}
	
	double acs()	{return acs_;}
	double acns()	{return acns_;}
	
	double ZCluster() {return Z;}
	double Xratlers() {return Xratlers_;}
	
	//double Xforce_simple() {return Xforce_simple_;}
	//double Xforce_double() {return Xforce_double_;}
	//double Xforce_triple() {return Xforce_triple_;}
	
	interdof* catchInterdof(unsigned int k)   { return linterdof_[k]; }
	
	double MeanRadial(unsigned int i) const {return MeanRadial_[i];}
	vector<double> MeanRadial() {return MeanRadial_;}
	
	double MeanOrtho(unsigned int j) const {return MeanOrtho_[j];}
	vector<double> MeanOrtho() {return MeanOrtho_;}
	

	vector<double> connectivity();
	vector<double> rankStat();
	vector<double> contactTypes();
	vector<double> contactDouble();
	void contact_connectivity();
	
	//void PDFFN_Cluster();
	//void PDFFT_Cluster();
	void Ptheta_Cluster(unsigned int,unsigned int, unsigned int);
	void contact_Stat();
	
	void MeanForce();

	void Stress_Cluster(Probe &);
	void Stress_ClusterSimple(Probe &);
	void Stress_ClusterNsimple(Probe &);
	
	//void Fabric_Cluster(Probe &); fabrique de texture
	void frAnisoInProbe(Probe &);
	void simpleContactFabric(Probe &);
	void NonsimpleContactFabric(Probe & prb);
	void ContactFabric(Probe & prb);
	void ContactForceRatio();
};

	

#endif //_clusterA_hpp
