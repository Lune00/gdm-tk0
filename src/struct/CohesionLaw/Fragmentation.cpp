#ifndef _Fragmentation_h
#define _Fragmentation_h

#include "cohesionLaw.hpp"


class Fragmentation : public cohesionLaw
{
protected: 
	unsigned int nsite_;
	double force_;
	vector< vector<unsigned int> > ActiveCohesion_;
public:
	
	double dAct( inter2d * inter ) {return 0.001*inter->first()->sizeVerlet();}
	
	void read( istream & in )   
	{
		in >> force_;
		cout<<"Force cohésion "<<force_<<endl;
		string filename;
		in >> filename;
 		read_data(filename.c_str());
	}

	void read_data(const char* name)
	{
		ifstream datafile(name);
  		if(!datafile)
		{
    		cerr << "@Simulation::read_data, cannot open file " << name << endl;
    		return;
   		}
    	
		string type;
    	unsigned int Npolyg;
  		vector<unsigned int> Active;
    
   		datafile >> type;
		while(datafile)
    	{
    		if (type == "Npolyg") datafile >> Npolyg;
			if(type=="nsite") datafile >> nsite_;
   			if (type=="Voisin")
   			{
   				Active.clear();
   		 		unsigned int voisin, nb;
   		 		datafile >> nb;
   		 		for(unsigned int i=0; i<nb; i++) {datafile >> voisin; Active.push_back(voisin);}
   		 		ActiveCohesion_.push_back(Active);
   			}
			datafile >> type;
		}
    	
   		cout<<"Npolyg:="<<Npolyg<<"   nsite:="<<nsite_<<endl;	
	}

	double fco( inter2d * inter )
	{
		unsigned int first, second;
		double force=0;
		first=inter->first()->id();
		second=inter->second()->id();
		
		for(unsigned int i=0; i<ActiveCohesion_[first].size(); i++)
		if(ActiveCohesion_[first][i]==second) force=force_;
		return force;
	} 
	
	vector< vector<unsigned int> > ActiveCohesion() {return ActiveCohesion_;}
	
	unsigned int nsite() {return nsite_;}
	
	double resistance() const {return force_;}
	
	Fragmentation():force_(0){ name_="fragmentation";}	
	
	~Fragmentation() {}
	
};

#endif // _Fragmentation_h


