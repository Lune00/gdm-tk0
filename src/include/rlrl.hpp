#ifndef _rlrl_hpp
#define _rlrl_hpp

#include "body2d.hpp"
#include "rline.hpp"
#include "inter2d.hpp"

//! \brief Interaction between two rlines
//! \author V. Richefeu
class rlrl : public inter2d
{
  
private:
  
  rline* i_;
  rline* j_;
  
  double facn0_, facn1_, facn2_, facn3_, facn4_;
  double fact0_, fact1_, fact2_, fact3_, fact4_;
  
  rlrl(rline* i,rline* j) : inter2d(), i_(i), j_(j) 
    {
    facn0_ = facn1_ = facn2_ = facn3_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    }
  
  rlrl() { facn0_ = facn1_ = facn2_ = facn3_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
  }
  rlrl(inter2d &) { }
  
  // only inter2d::factory can create an inter2d child
  friend class inter2d;	
	
public:
		
  ~rlrl() {}
	
  rline* first()  { return i_; }
  rline* second() { return j_; }
  
  body2d* lexifirst()  { return i_; }
  body2d* lexisecond() { return j_; }
  
  void findProxApex();
	
  double Dist();
  double Overlap();
  void   Vel();
  void   Frame();
  void   Pos();
  void   Kin();
  void   Res();
  void   Res(const double, const double);
  void   CDcoeff();
  double An(const double dt);
  double At(const double dt);
  void   read(istream & is);
  void   write(ostream & os);
  void   writeMGP(ostream & os);
double fx()   { return (fn()*nx_ + ft()*tx_);}
double fy()   { return (fn()*ny_ + ft()*ty_);}
	
};

#endif // _rlrl_hpp
