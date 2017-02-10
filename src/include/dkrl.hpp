#ifndef _dkrl_hpp
#define _dkrl_hpp

#include "body2d.hpp"
#include "disk.hpp"
#include "rline.hpp"
#include "inter2d.hpp"

//! \brief Interaction between a disk and a rline
//! \author V. Richefeu
class dkrl : public inter2d
{
	
  disk*  i_;
  rline* j_;
  
  body2d* lexi0_;
  body2d* lexi1_;

double overlap_;

  // These variables are used for stocking some predetermined factors
  double facn0_, facn1_, facn2_, facn3_, facn4_;
  double fact0_, fact1_, fact2_, fact3_, fact4_;
  double facs0_, facs1_, facs2_;

  dkrl() {  facn0_ = facn1_ = facn2_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
	type()=_type_dkrl;
}
  dkrl(inter2d &) {  facn0_ = facn1_ = facn2_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
	type()=_type_dkrl;
}
  dkrl(disk* i,rline* j, bool lexiswap) : inter2d(), i_(i), j_(j) 
    { 
    facn0_ = facn1_ = facn2_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
	type()=_type_dkrl;
    
    // voir dkdkP...
    if (lexiswap) { lexi0_ = j; lexi1_ = i; }
    else          { lexi0_ = i; lexi1_ = j; }
    }
  
  friend class inter2d;

public:
	
  ~dkrl() {}
	
  disk*  first()  { return i_; }
  rline* second() { return j_; }
  
  void plug(body2d* b1, body2d* b2) 
    {
    i_ = dynamic_cast<disk*>(b1);
    j_ = dynamic_cast<rline*>(b2);

    if(i_->id() < j_->id())
      {
      lexi0_ = b1;
      lexi1_ = b2;
      }
    else
      {
      lexi0_ = b2;
      lexi1_ = b1;
      }
    }

  body2d* lexifirst()  { return lexi0_; }
  body2d* lexisecond() { return lexi1_; }
  
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
  
  void clearForceAndMoment() { fn_ = ft_ = frot_ = 0.0;vn_ = vt_ = vrot_ = 0.0; }
  
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

  double   Vbranchx() const {return (i_->x() - j_->x());}
  double   Vbranchy() const{return (i_->y() - j_->y());}
double fx()  const { return (fn()*nx_ + ft()*tx_);}
double fy()  const { return (fn()*ny_ + ft()*ty_);}
};

#endif // _dkrl_hpp
