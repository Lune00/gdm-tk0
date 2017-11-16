#include "sample.hpp"

using namespace std;

vector<dof*> Sample::read(istream & is)
{
	string type;
	int n=0;
	vector<dof*> ldof;
	is >> type;	
	while(is)
	{
		if(type == "file")
		{
			string filename;
			is >> filename;
			ifstream inc(filename.c_str());
			if(inc)
			{
				Sample::read(inc);
			}
			else
			{
				cerr << "@Sample::read, cannot open file " << filename << endl;
			}
		}

		if(type == "}") break;	

		if(type == "Cluster{")
		{
			ldof.push_back(new dof);
			ldof.back()->id()=ldof.size()-1;
			is>>type;
			while( is)
			{
				if(type == "}") break;
				else if (type == "area")
					is>>ldof.back()->Area();
				else
				{
				lbody_.push_back(body2d::factory(type));
				lbody_[n]->read(is);
				lbody_[n]->id() = n;
				
				ldof.back()->plugBody(lbody_[n]);
				++n;
				}
				is>>type;
			}
		}
		else
		{
			lbody_.push_back(body2d::factory(type));
			lbody_[n]->read(is);
			lbody_[n]->id() = n;
			++n;
		}

		is >> type;
	}
	cout<<"Lecture termine : "<<lbody_.size()<<endl;
	return (ldof);
}

void Sample::write(const char * fname)
{
	ofstream fsample(fname);

	for (unsigned int i=0;i<lbody_.size();++i)
	{
		if(lbody_[i]->bodyDof() != NULL )
		{
			fsample<<"Cluster{"<<endl;
			i+=lbody_[i]->bodyDof()->writeSpl(fsample);
			if( lbody_[i]->bodyDof()->Area()!= -1)
				fsample<<"area "<<lbody_[i]->bodyDof()->Area()<<endl;
			fsample<<"}"<<endl;
		}
		else
			lbody_[i]->write(fsample);
			//this->body(i)->write(fsample);
	}
}


void Sample::write(ostream & os)
{	
	for (unsigned int i=0;i<lbody_.size();++i)
	{
		if(lbody_[i]->bodyDof() != NULL )
		{
			os<<"Cluster{"<<endl;
			i+=lbody_[i]->bodyDof()->writeSpl(os);			
			os<<"}"<<endl;
		}
		else
		{
		//	if(i<2) cout<<lbody_[0]->grp()<<"	Dien tich "<<i<<"	"<<lbody_[i]->Area()<<endl;
			lbody_[i]->write(os);
			//this->body(i)->write(os);
		}
			
	}	
}

void Sample::removeBody(body2d * rm)
{

	vector< body2d*>::iterator itspl;
		for(itspl = this->lbody().begin() ; itspl != this->lbody().end() ; )
		{
			if( (*itspl)==rm)
				itspl=this->lbody().erase( itspl);
			else
				itspl++;
		}
	}

//Assigne la masse a chaque particule à partir de la densité.
void Sample::fill(double density)
{
	for (unsigned int i=0;i<lbody_.size();++i)
		lbody_[i]->Fill(density);
}


//Calcule xmin,xmax,ymin,ymax de l'échantillon
void Sample::updateBoundaries()
{
	unsigned int N = lbody_.size();
	if(N == 0) 
	{
		xmin_ = xmax_ = ymin_ = ymax_ = 0.0;
		return;
	}

	double val;
    //cout<<"update boundaries()"<<endl;
	xmin_ = lbody_[0]->x();
	xmax_ = lbody_[0]->x();
	ymin_ = lbody_[0]->ymin();
	ymax_ = lbody_[0]->ymax();

	for (unsigned int i=1 ; i<N ; ++i)
	{
		xmin_ = xmin_ > (val = lbody_[i]->xmin()) ? val : xmin_;
		xmax_ = xmax_ < (val = lbody_[i]->xmin()) ? val : xmax_;
		ymin_ = ymin_ > (val = lbody_[i]->ymin()) ? val : ymin_;
		ymax_ = ymax_ < (val = lbody_[i]->ymax()) ? val : ymax_;
	}
}


// ajouter calcul rmoy_
void Sample::radiusExtrema( unsigned int start) //updateRadii()
{
	rmoy_=0.;
	unsigned int N = lbody_.size();
	if(N == 0) 
	{
		rmin_ = rmax_ = 0.0;
		return;
	}

	double val=0;

	rmin_ = rmax_ =lbody_[start]->sizeVerlet(); // Rayon du body 4


	for (unsigned int i=start+1 ; i<N ; ++i)
	{
		rmin_ = rmin_ > (val = lbody_[i]->sizeVerlet()) ? val : rmin_;
		rmax_ = rmax_ < (val = lbody_[i]->sizeVerlet()) ? val : rmax_;
		rmoy_+=lbody_[i]->sizeVerlet();
	} 
	rmoy_/=(double) (N);
}



void Sample::definePeriodicity(double bandwidth)
{
    
    if( bandwidth < 1. ) bandwidth = 2 * rmax_;
   
		
		isMonoPeriodic_ = true;
		leftBoundary_   = xmin_;
		rightBoundary_  = xmax_;
		boundWidth_     = (rightBoundary() - leftBoundary());
		bandWidth_      = bandwidth * rmax_ ;


}

//Remet a jour les bandes périodiques
void Sample::updateBands()
{
	unsigned int N = lbody_.size();
	double      lb = leftBoundary();
	double      rb = rightBoundary();
	double      bw = bandWidth(); 

	leftband().clear();
	rightband().clear();

	for (unsigned int i=0;i<N;++i)
	{
		if (lbody_[i]->x() <= (lb + bw))
		{
			leftband_.push_back(i);
		}

		if (lbody_[i]->x() >= (rb - bw))
		{
			rightband_.push_back(i);
		}
	}
}

void Sample::translate()
{
	unsigned int NR = rightband_.size();
	unsigned int NL = leftband_.size();
	double lb  = leftBoundary();
	double rb  = rightBoundary();
	double bow = boundWidth();
	unsigned int taille = 1000;
	vector <dof*> ldof(taille);//WARNING NUMBER OF DOF IS HERE LIMITED.......................................
	unsigned int i;
	for ( i=0; i< taille; ++i)
		ldof[i]=NULL;
		
	for (i = 0 ; i < NR ; ++i)
	{
		if (body( rightband(i))->bodyDof()==NULL)
		{
			if (body( rightband(i))->x() > rb) 
				body(rightband(i))->x() -= bow;
		}
		else 
		{
			ldof[body( rightband(i))->bodyDof()->id() ]=body( rightband(i))->bodyDof();
		}
	}
	for (i = 0 ; i < NL ; ++i)
	{
		if (body(leftband(i))->bodyDof()==NULL)
		{
			if (body(leftband(i))->x() < lb) 
				body(leftband(i))->x() += bow;
		}
		else
		{
			ldof[body( leftband(i))->bodyDof()->id() ]=body( leftband(i))->bodyDof();
		}
	}

	
	for ( i=0; i< taille; ++i)
	{
		if ( ldof[i]==NULL ) continue;

		if (  ldof[i]->isPeriodic() )
		{
			ldof[i]->translateBodies( rb, lb, bow );
			//notTranslate=false;
		}
		else
		{
			if ( ldof[i]->leftBody()->xmin() > rb) 
				{
					cout<<" translation negative dof "<<i<<endl;
					ldof[i]->translate(-bow) ;
				}
			else if ( ldof[i]->rightBody()->xmax() < lb) 
				{
					cout<<" translation positive dof "<<i<<endl;
					ldof[i]->translate(bow) ;
				}
		}
	}


}

void Sample::SortByHeight() // from bottom to top
{
	std::stable_sort(lbody_.begin(), lbody_.end(), body2d::compareHeight() );
}	

void Sample::swapBodies(unsigned int a,unsigned int b)
{
	body2d *tmp;
	tmp = body(a);
	lbody_[a] = body(b);
	lbody_[b] = tmp;
}
