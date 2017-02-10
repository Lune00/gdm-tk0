#ifndef _connectivity_A_hpp
#define _connectivity_A_hpp

#include "sample.hpp"
#include "network.hpp"
#include "probe.hpp"



using namespace std;

class connectivity_A
{
private : 
	
	//donner un vecteur de inter2D* en argument et une liste de fonction plutot qu'une classe ???
 Probe * prb_;
 Sample * spl_;
 Network * nwk_;
 
	vector< inter2d * > contacts_;// Interactions in contact
	vector< bool > activatedBodies_;// Bodies with non zero applied forces
	
	//All measures are done in the probe...
	vector< unsigned int > Ncloc_;//Number of contacts per activated bodies
	vector< unsigned int > Npi_;//Number of activated particles with i contacts (normally begin to 2)
	double z_;//Global coordination number
	unsigned int Noc_;//total number of contacts
	unsigned int N;
	
	
	unsigned int Nosc_;//total number of sliding contacts
	unsigned int Norc_;//total number of rolling contacts
	
	
	double fracRa_;//fraction of rattlers bodies

public :
	
	connectivity_A(){};
	
	connectivity_A(Probe * prb,Sample * spl,Network * nwk):
	prb_(prb),spl_(spl),nwk_(nwk) {	};
	
	~connectivity_A(){};

	void plug (Probe * prb,Sample * spl,Network * nwk){ prb_=prb;spl_=spl;nwk_=nwk;}
	
	void check();//check for contacts and linked bodies in probe
	
	double Z();//evaluation of z_
	
//	pointSet * sizeCorrelation( );//correlation between size of particles and number of contacts
	
	
	double   z() const {return z_;}
	double & z()       {return z_;}
	
	vector < double > Pq();//fraction of particles with q contacts 
	
	unsigned int   Npi( unsigned int i ) const 
	{ 
		if (i< Npi_.size() ) return (Npi_[i]);
		else return (0);
	}
	
	unsigned int   Noc() const {return Noc_;}
	unsigned int & Noc()       {return Noc_;}
	
	unsigned int   Nosc() const {return Nosc_;}
	unsigned int & Nosc()       {return Nosc_;}
	

};

#endif // _connectivity_A_hpp
