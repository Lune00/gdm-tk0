#ifndef _dkdkP_h
#define _dkdkP_h

#include "body2d.hpp"
#include "disk.hpp"
#include "inter2d.hpp"

//! \brief Interaction between a disk and an periodic-image of disk
//! \author C. Voivret
class dkdkP : public inter2d
{
	
  disk* i_;
  disk* j_;
    
  double P_;

  body2d* lexi0_;
  body2d* lexi1_;
  
  double facn0_,facn1_,facn2_,facn3_,facn4_;
  double fact0_,fact1_,fact2_,fact3_,fact4_;
  double facs0_,facs1_,facs2_;

  dkdkP() { type_=_type_dkdkP;}
  dkdkP(inter2d &) { type_=_type_dkdkP; }
  dkdkP(disk* i,disk* j,double P) : inter2d(), i_(i), j_(j), P_(P) 
  {
    facn0_ = facn1_ = facn2_ = facn3_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
    
    lexiOrder();
	type_=_type_dkdkP;
  }
  
  friend class inter2d;  

public:
		
  ~dkdkP() { }
	
  disk* first()  { return i_; }
  disk* second() { return j_; }
  
  void plug(body2d* b1, body2d* b2) 
    {
    i_ = dynamic_cast<disk*>(b1);
    j_ = dynamic_cast<disk*>(b2);
	lexiOrder();
    }

  body2d* lexifirst()  { return lexi0_; }
  body2d* lexisecond() { return lexi1_; }

double & P()          { return P_; }
  double   P()    const { return P_; }
  void lexiOrder()
  {
  if (i_->id() < j_->id())
      { lexi0_ = i_; lexi1_ = j_;}
    else
      { lexi0_ = j_; lexi1_ = i_;}
	}
	
  double Dist();
  double Overlap();
  void   Frame();
  void   Frame2();
  bool   Activate();
  void   Kin();
  void   Res();
  void   Res(const double dfn, const double dft, const double dfs = 0.0);
  void   CDcoeff();
  double An(const double dt);
  double At(const double dt);
  double As(const double dt);
  void   read(istream & is, unsigned int * Id1, unsigned int * Id2);
  void   write(ostream & os);
  void   writeMGP(ostream & os);	
  
  void clearForceAndMoment() { fn_ = ft_ = frot_ = 0.0; }
  
  double & x()        { return x_;  }
  double   x()  const { return x_;  }
  double & y()        { return y_;  }
  double   y()  const { return y_;  }	
  
  double & vn()       { return vn_; }
  double   vn() const { return vn_; }
  double & vt()       { return vt_; }
  double   vt() const { return vt_; }
  
  double & fn()       { return fn_; }
  double   fn() const { return fn_; }
  double & ft()       { return ft_; }
  double   ft() const { return ft_; }

  double   Vbranchx() const{return (i_->x() - j_->x() - P_ );}
  double   Vbranchy() const{return (i_->y() - j_->y());}

double fx()  const { return (fn()*nx_ + ft()*tx_);}
double fy()  const { return (fn()*ny_ + ft()*ty_);}
};

#endif // _dkdkP_h
