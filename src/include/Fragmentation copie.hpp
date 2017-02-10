#ifndef _Fragmentation_hpp
#define _Fragmentation_hpp

#include "cohesionLaw.hpp"
#define MIN(A,B) ((A)<(B) ? (A):(B))

class Fragmentation : public cohesionLaw
{
protected: 
	unsigned int nsite_;
	double force_;
	vector< vector<unsigned int> > ActiveCohesion_;
public:
	
	double dAct( inter2d * inter ) {return 0;}
	//{return 0.000001*MIN(inter->first()->sizeVerlet(),inter->second()->sizeVerlet());}
	
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
		 for(unsigned int i=0;i<ActiveCohesion_.size();i++)
		{
			for(unsigned int j=0;j<ActiveCohesion_[i].size();j++) cout<<ActiveCohesion_[i][j]<<"|";
			cout<<endl;
    	}
	}

	double fco( inter2d * inter )
	{
		unsigned int first, second;
		double force=0;
		first=inter->first()->id()-4;
		second=inter->second()->id()-4;
		
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

#endif // _Fragmentation_hpp


