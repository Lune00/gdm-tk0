#ifndef _rectangularProbe_hpp
#define _rectangularProbe_hpp

#include "probe.hpp"

//! \brief Rectangular Probe
//! \author V. Richefeu
class rectangularProbe : public Probe
{	
protected:
  
double xc_,yc_;
double hh_,hl_;//half height, half length
  
public:
  
  rectangularProbe(): xc_(0.0), yc_(0.0), hh_(1.0),hl_(1.0) { }
  rectangularProbe(double x, double y, double hh, double hl): xc_(x), yc_(y), hh_(hh), hl_(hl) { }
  ~rectangularProbe() { }
  
  
  
  double area() {return (4*hh_*hl_);}
  
  bool containDofCenter(dof * i)
  {
    double dx = fabs(i->mcx() - xc_);
    double dy = fabs(i->mcy() - yc_);
    return ( dx < hl_ && dy < hh_ );

  }
  
  bool containCenter (body2d * b)
  {
    double dx = fabs(b->x() - xc_);
    double dy = fabs(b->y() - yc_);
    return ( dx < hl_ && dy < hh_ );
  }
  
  bool contain (inter2d * k)
  {
    double dx = fabs(k->x() - xc_);
    double dy = fabs(k->y() - yc_);
    return ( dx < hl_ && dy < hh_ );
  }
  
  /*
  bool intersection ( body2d * b)
  {
	if( xc_ + hl_ >= b->xmin()) return true;
	if( xc_ - hl_ <= b->xmax()) return true;
	
	if( yc_ + hh_ >= b->ymin()) return true;
	if( yc_ - hh_ <= b->ymax()) return true;
	
	return false;	
  
  }*/
  
  // Corriger par DH NGUYEN
  bool intersection ( body2d * b)
  {
	if( xc_ + hl_ <= b->xmin()) return false;
	if( xc_ - hl_ >= b->xmax()) return false;
	
	if( yc_ + hh_ <= b->ymin()) return false;
	if( yc_ - hh_ >= b->ymax()) return false;
	
	return true;	
  
  }
  
  
  bool containEntireBody ( body2d * b )
  {
	if(
	    xc_ + hl_ >= b->xmax()  
	 && xc_ - hl_ <= b->xmin() 
	 && yc_ + hh_ >= b->ymax()
	 && yc_ - hh_ <= b->ymin() 
		) return true;
				
	return false;
	
  }

bool probeInsideBody ( body2d * b )
  {
	cout<<" @ rectangularProbe : probeInsideBody not defined "<<b->id()<<endl;
	return false;
  }

  double & x()          { return xc_; }
  double   x()    const { return xc_; }
  double & y()          { return yc_; }
  double   y()    const { return yc_; }
  double & hh()          { return hh_; }
  double   hh()    const { return hh_; } 
  double & hl()          { return hl_; }
  double   hl()    const { return hl_; }
};

#endif // _circularProbe_hpp
