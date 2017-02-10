#ifndef _Fragmentation_hpp
#define _Fragmentation_hpp

#define MIN(A,B) ((A)<(B) ? (A):(B))

class Fragmentation
{
protected: 
	unsigned int Npolyg;
	unsigned int nsite_;
	unsigned int ngrain;
	long int numfissure_;
	long int ncontact_;
	double frict_; //coefficient de friction
	double cohesion_; //pression cohésion`
	double epsilon_;
	vector< vector<unsigned int> > ActiveCohesion_;
	vector <bool> lgrain_;//rentrer état de particule fracture ou non //tra ve gia tri dung neu hat i ban dau da co vet nut
	vector <bool> casse_; //rentrer true si cette cellule de Voronoi est cassé //tra ve gia tri dung neu Vovonoi i thuoc particule da casse
	vector< vector<unsigned int> > Cluster_;
public:

	Fragmentation() {cohesion_=0;frict_=0;epsilon_=0.001;nsite_=0;numfissure_=0;ncontact_=0;}	
	
	~Fragmentation() {}
	
	void read( istream & is )   
	{
		is >> cohesion_ >> frict_;
		string filename;
		is >> filename;
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
    	unsigned int Ncluster;
  		vector<unsigned int> Active;
  		vector<unsigned int> ListPolygone;
    
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
    	
    	ngrain=Npolyg/nsite_;
    	for (unsigned int i=0; i<ngrain; i++) lgrain_.push_back(false);
    	for (unsigned int i=0; i<ngrain; i++) casse_.push_back(false);
    	
   		cout<<"Npolyg:="<<Npolyg<<"   nsite:="<<nsite_<<endl;
		for(unsigned int i=0;i<ActiveCohesion_.size();i++)
		{
			ncontact_+=ActiveCohesion_[i].size();
			for(unsigned int j=0;j<ActiveCohesion_[i].size();j++) cout<<ActiveCohesion_[i][j]<<"|";
			cout<<endl;
    	}
    	ncontact_/=2;
    	Ncluster=Npolyg/nsite_;
    	for(unsigned int i=0; i<Ncluster; i++)
    	{
    		ListPolygone.clear();
    		for(unsigned int j=0; j<nsite_;j++) ListPolygone.push_back(i*nsite_+j);
    		Cluster_.push_back(ListPolygone);
    	}	
	}
	
	void write(ofstream & fra)
	{
		fra<<"Fragmentation{"<<endl;
		fra<<"Cohesion "<<cohesion_<<" Frict "<<frict_<<endl;
		fra<<"Npolyg "<<Npolyg<<" Nsite "<<nsite_<<endl;
		fra<<"Nfissure "<<numfissure_<<endl;
		
		for(unsigned int i=0; i<ActiveCohesion_.size(); i++)
		{
			fra<<"Voisin "<<ActiveCohesion_[i].size()<<" ";
			for(unsigned int j=0; j<ActiveCohesion_[i].size(); j++) fra<<ActiveCohesion_[i][j]<<" ";
			fra<<endl;
		}
		
		for(unsigned int i=0; i<Cluster_.size(); i++)
		{
			fra<<"FCluster "<<Cluster_[i].size()<<" ";
			for(unsigned int j=0; j<Cluster_[i].size(); j++) fra<<Cluster_[i][j]<<" ";
			fra<<endl;
		}
		
		fra<<"Casse"<<endl;
		for(unsigned int i=0; i<ngrain; i++)
		{
			fra<<casse_[i]<<" ";
			if (i%50==49) fra <<endl;
		}
		fra<<endl;
		
		fra<<"Fracture"<<endl;
		for(unsigned int i=0; i<ngrain; i++)
		{
			fra<<lgrain_[i]<<" ";
			if (i%50==49) fra <<endl;
		}
		fra<<endl;
		
		fra<<"}"<<endl<<endl;
	}
	
	void load_history(istream & fra)
	{
		string token;
		bool etat;
		vector<unsigned int> Active;
  		vector<unsigned int> ListPolygone;
		
		fra >> token;
		while(fra)
		{
			if (token == "Cohesion") fra >> cohesion_;
			else if (token == "Frict") fra >> frict_;
			else if (token == "Npolyg") fra >> Npolyg;
			else if (token == "Nsite") fra >> nsite_;
			else if (token == "Nfissure") fra >> numfissure_;
			else if (token == "Voisin")
			{
				Active.clear();
   		 		unsigned int voisin, nb;
   		 		fra >> nb;
   		 		for(unsigned int i=0; i<nb; i++) {fra >> voisin; Active.push_back(voisin);}
   		 		ActiveCohesion_.push_back(Active);
			}
			else if (token == "FCluster") 
			{
				ListPolygone.clear();
				unsigned int voisin, nb;
   		 		fra >> nb;
   		 		for(unsigned int i=0; i<nb; i++) {fra >> voisin; ListPolygone.push_back(voisin);}
   		 		ActiveCohesion_.push_back(Active);
				Cluster_.push_back(ListPolygone);
			}
			else if (token == "Casse")
			{
				casse_.clear();
				for(unsigned int i=0; i<ngrain; i++) { fra >> etat; casse_.push_back(etat);}
			}
			else if (token == "Fracture")
			{
				lgrain_.clear();
				for(unsigned int i=0; i<ngrain; i++) { fra >> etat; lgrain_.push_back(etat);}
			}
			else if (token == "}") break;
		
			fra >> token;
		}		
	}


	double fco(inter2d * inter )
	{
		unsigned int first, second;
		double fco=0;
		first=inter->first()->id()-4;
		second=inter->second()->id()-4;
		
		unsigned int i=0;
		bool logic=true;
		
		while (logic && (i<ActiveCohesion_[first].size()))
		{
			if(ActiveCohesion_[first][i]==second)
			{
				fco=0.5*cohesion_*inter->longeur();
				logic=false;
			}
			else i++;
		}
		return fco;
	} 
	
	double frict(inter2d * inter )
	{
		unsigned int first, second;
		double frict=0;
		first=inter->first()->id()-4;
		second=inter->second()->id()-4;
		
		unsigned int i=0;
		bool logic=true;
		
		while (logic && (i<ActiveCohesion_[first].size()))
		{
			if(ActiveCohesion_[first][i]==second)
			{
				frict=frict_;
				logic=false;
			}
			else i++;
		}
	
		return frict;
	}

	vector< vector<unsigned int> > & ActiveCohesion()       {return ActiveCohesion_;}
	vector< vector<unsigned int> >   ActiveCohesion() const {return ActiveCohesion_;}
	
	vector< vector<unsigned int> > & Cluster()       {return Cluster_;}
	vector< vector<unsigned int> >   Cluster() const {return Cluster_;}
	
	vector< bool> & Casse()       {return casse_;}
	vector< bool >   Casse() const {return casse_;}
	
	vector <bool> & lgrain() 		  {return lgrain_;}
	vector <bool>   lgrain()	const {return lgrain_;}
	
	unsigned int nsite() const {return nsite_;}
	long int & numfissure() {return numfissure_;}
	long int numfissure() const {return numfissure_;}
	long int & ncontact() {return ncontact_;}	
	
	double sigmac() {return cohesion_;}
	double epsilonc() {return epsilon_; }
	double frict() {return frict_; } 
	
	double res_traction(inter2d * inter)  {return -0.5*cohesion_*inter->longeur();}
	
	double res_cisaillement(inter2d * inter) {return 0.5*frict_*cohesion_*inter->longeur();}
	
	bool fracture(inter2d * inter)
	{
		double critere1,critere2;
		
		if(inter->rang()>1)
		{
			inter->current()=0;
			critere1=inter->fn()/this->res_traction(inter)+fabs(inter->ft())/this->res_cisaillement(inter);
		
			inter->current()=1;
			critere2=inter->fn()/this->res_traction(inter)+fabs(inter->ft())/this->res_cisaillement(inter);
			
			inter->current()=0;
			if ((critere1>1.)&&(critere2>1.)) return true ;
			else return false;
		}
		else return true;		
	}
};

#endif // _Fragmentation_hpp


