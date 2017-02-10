#ifndef _pgrl_hpp
#define _pgrl_hpp

#include "body2d.hpp"
#include "polyg.hpp"
#include "rline.hpp"
#include "inter2d.hpp"

//! \brief Interaction between a disk and a rline
//! \author V. Richefeu
class pgrl : public inter2d
{
	
  polyg*  i_;
  rline* j_;
  
  body2d* lexi0_;
  body2d* lexi1_;
double overlap_;
  // These variables are used for stocking some predetermined factors
  double facn0_, facn1_, facn2_, facn3_, facn4_;
  double fact0_, fact1_, fact2_, fact3_, fact4_;
  double facs0_, facs1_, facs2_;
  double facn0_2_, facn1_2_, facn2_2_, facn3_2_, facn4_2_;
  double fact0_2_, fact1_2_, fact2_2_, fact3_2_, fact4_2_;


  double x2_ , y2_;   // Position of the second contact point
  double fn2_, ft2_;
  double vn2_, vt2_;

  pgrl() {facn0_ = facn1_ = facn2_ = facn3_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
    facn0_2_ = facn1_2_ = facn2_2_ = facn3_2_ = facn4_2_ = 0.0;
    fact0_2_ = fact1_2_ = fact2_2_ = fact3_2_ = fact4_2_ = 0.0;

	type()=_type_pgrl; }
  pgrl(inter2d &) { facn0_ = facn1_ = facn2_ = facn3_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
    facn0_2_ = facn1_2_ = facn2_2_ = facn3_2_ = facn4_2_ = 0.0;
    fact0_2_ = fact1_2_ = fact2_2_ = fact3_2_ = fact4_2_ = 0.0;

	type()=_type_pgrl;
 }
  pgrl(polyg* i,rline* j, bool lexiswap) : inter2d(), i_(i), j_(j) 
    { 
    facn0_ = facn1_ = facn2_ = facn3_ = facn4_ = 0.0;
    fact0_ = fact1_ = fact2_ = fact3_ = fact4_ = 0.0;
    facs0_ = facs1_ = facs2_ = 0.0;
    facn0_2_ = facn1_2_ = facn2_2_ = facn3_2_ = facn4_2_ = 0.0;
    fact0_2_ = fact1_2_ = fact2_2_ = fact3_2_ = fact4_2_ = 0.0;

	type()=_type_pgrl;
    
    // voir dkdkP...
    if (lexiswap) { lexi0_ = j; lexi1_ = i; }
    else          { lexi0_ = i; lexi1_ = j; }
    }
  
  friend class inter2d;

public:
	
  ~pgrl() {}
	
  polyg*  first()  { return i_; }
  rline* second() { return j_; }
  
  void plug(body2d* b1, body2d* b2) 
    {
    i_ = dynamic_cast<polyg*>(b1);
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
double Overlap(){return overlap_;}
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
  
  void clearForceAndMoment() 
	{ 
	fn_  = ft_  = frot_ = 0.0;
    fn2_ = ft2_ = 0.0;
 	}
  
  double & x()
    {
    if ( rank_ == 2 && current_ == 1) return x2_;
    return x_;
    }
  
  double  x()  const 
    {
      if (rank_ == 2 && current_ == 1) return x2_;
      return x_;
    }
  
  double & y() 
    {
    if ( rank_ == 2 && current_ == 1) return y2_;
    return y_;
    }
  double  y()  const 
    {
      if (rank_ == 2 && current_ == 1) return y2_;
      return y_;
    }
  
  double & vn() 
    {
    if (rank_ == 2 && current_ == 1) return vn2_;
    return vn_;
    }
  
  double vn() const
    {
      if (rank_ == 2 && current_ == 1) return vn2_;
      return vn_;
    }
  
  double & vt() 
    {
    if (rank_ == 2 && current_ == 1) return vt2_;
    return vt_;
    }
  
  double vt() const 
    {
      if (rank_ == 2 && current_ == 1) return vt2_;
      return vt_;
    }

  double & fn() 
    {
    if (rank_ == 2 && current_ == 1) return fn2_;
    return fn_;
    }
  
  double fn() const
    {
      if (rank_ == 2 && current_ == 1) return fn2_;
      return fn_;
    }
  
  double & ft()
    {
    if (rank_ == 2 && current_ == 1) return ft2_;
    return ft_;
    }
  
  double ft() const
    {
      if (rank_ == 2 && current_ == 1) return ft2_;
      return ft_;
    }

	double   Vbranchx() const{return (i_->x() - j_->x());}
	double   Vbranchy() const{return (i_->y() - j_->y());}
	
	double fn1() const
	{
		return fn_;
	}
	
	double fn2() const
	{
		if( rank_==2) return fn2_;
		else return 0.;
	}

	double ft1() const
	{
		return ft_;
	}
	double ft2() const
	{
		if (rank_==2) return ft2_;
		else return 0.;
	}
	double fx() const   
	{
		double fx=0;
		fx = fn1()*nx_ + ft1()*tx_;
		if (rank_==2)
		{
		fx+=fn2()*nx_ + ft2()*tx_;
		}
		return (fx);
	}

	double fy() const  
	{  
		double fy=0;
		fy = fn1()*ny_ + ft1()*ty_;
		if (rank_==2)
		{
		fy+=fn2()*ny_ + ft2()*ty_;
		}
		return (fy);
	}
};

#endif // _pgrl_hpp
