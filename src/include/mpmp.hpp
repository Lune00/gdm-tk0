#ifndef _mpmp_hpp
#define _mpmp_hpp

#include "body2d.hpp"
#include "massPoint.hpp"
#include "inter2d.hpp"

//! \brief Interaction between two massPoints
//! \author V. Richefeu
class mpmp : public inter2d
{
	
  massPoint* i_;
  massPoint* j_;

  mpmp() { }
  mpmp(inter2d &) { }
  mpmp(massPoint* i, massPoint* j) : inter2d(), i_(i), j_(j) 
  {
  }
  
  friend class inter2d;  

	public:
		
  ~mpmp() { }
	
  massPoint* first()  { return i_; }
  massPoint* second() { return j_; }
  
  void plug(body2d* b1, body2d* b2) 
    {
    i_ = dynamic_cast<massPoint*>(b1);
    j_ = dynamic_cast<massPoint*>(b2);
    }

  body2d* lexifirst()  { return i_; }
  body2d* lexisecond() { return j_; }
	
  double Dist();
double Overlap(){return 0;}
  void   Frame();
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
 double fx()   { return (fn()*nx_ + ft()*tx_);}
 double fy()   { return (fn()*ny_ + ft()*ty_);}	
};

#endif // _mpmp_h
