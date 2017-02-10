#include "connectivity_A.hpp"

void connectivity_A::check( )
{
	contacts_.clear();
	activatedBodies_.clear();
	Ncloc_.clear();
	
	activatedBodies_.resize( spl_->lbody().size(),false);
	
	Ncloc_.resize( spl_->lbody().size(),0);
	
	cout<<activatedBodies_.size()<<" "<<Ncloc_.size()<<endl;
	unsigned int id1,id2;
	for(unsigned int i=0;i< nwk_->linter().size(); ++i )
	{
		if( prb_->contain( nwk_->inter(i) ))
		{
			if( nwk_->inter(i)->rang()>0 )
			{	
				contacts_.push_back(nwk_->inter(i));
				
				id1 = contacts_.back()->first()->id();
				id2 = contacts_.back()->second()->id();
				
				if( prb_->containCenter(contacts_.back()->first()))
				activatedBodies_[contacts_.back()->first()->id()]=true;
				
				if( prb_->containCenter(contacts_.back()->second()))
				activatedBodies_[contacts_.back()->second()->id()]=true;
				
				Ncloc_[id1 ]++;
				Ncloc_[id2 ]++;
				
			}
			else
			{
				
			}
			
		}
	}
	
	Npi_.clear();
	Npi_.resize(100,0);
	for (unsigned int i=0; i< activatedBodies_.size();++i)
	{
		if (activatedBodies_[i] )
		{
			Npi_[ Ncloc_[i] ]++;
			if( Ncloc_[i]==2) cout<<i<<endl;
		}
		 
	}
	
}

double connectivity_A::Z()
{
	unsigned int Nb=0,Nc=0,Nr=0;
	
	for (unsigned int i=0; i< activatedBodies_.size();++i)
	{
		if (activatedBodies_[i] )
		{
			if( Ncloc_[i]>1 )
			{	
				Nb++;
				Nc+=Ncloc_[i];
			}
			else
			{
				Nr++;
			}
		}
	}
	
	z_=(double) (Nc)/(double) (Nb);
	
	fracRa_= (double) (Nr)/(double) (Nr+Nb);
	
	return (z_);
}





