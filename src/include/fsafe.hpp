#ifndef _fsafe_hpp
#define _fsafe_hpp

#include "inter2d.hpp"

// pour le moment ne fonctionne que pour les "mono-contacts"
// \brief Force to be stocked
// \author V. Richefeu
class fsafe
{
  body2d *i_, *j_;
  double fn_, ft_;
  double fn2_,ft2_;
  unsigned int type_;
  unsigned int rank_;
	
 public:

  fsafe(inter2d *k) :
    i_(k->lexifirst()), j_(k->lexisecond()), 
    fn_(k->fn())  , ft_(k->ft()), 
	fn2_(0) , ft2_(0),
	type_( k->type() ),rank_(k->rang())
    { 
		//cout<<type_<<" "<<i_->id()<<" "<<j_->id()<<endl;
		if (rank_==2)
		{
			
		//	cout<<" double stock"<<endl;
			k->current()=1;
			fn2_=k->fn();
			ft2_=k->ft();
			k->current()=0;
			
		}
		
	/*	if(k->lexifirst()->id()<100  ) //&& k->lexifirst()->id()>-1 
		//if (fn_ != 0. || fn2_ != 0.)
		{
			cout<<"rank= "<<rank_<<" type "<<type_<<" "; 
			cout<<"stock "<< k->lexifirst()->id() <<" et "<< k->lexisecond()->id()<<" ";
			cout<<" 	fn1= "<<fn_<<" ft1= "<<ft_<<endl;
			//cout<<" 	fn2= "<<fn2_<<" ft2= "<<ft2_<<endl;
		} 
	*/	
	}

  body2d* first()  { return i_; }
  body2d* second() { return j_; }
  
  double & fn()       { return fn_; } 
  double   fn() const { return fn_; }
  double & ft()       { return ft_; }
  double   ft() const { return ft_; }
  
  double & fn2()       { return fn2_; } 
  double   fn2() const { return fn2_; }
  double & ft2()       { return ft2_; }
  double   ft2() const { return ft2_; }

  unsigned int & type()       { return type_; }
  unsigned int   type() const { return type_; }

  unsigned int & rang()       { return rank_; }
  unsigned int   rang() const { return rank_; }
	
};

#endif // _fsafe_hpp
