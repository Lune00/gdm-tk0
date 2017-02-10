#ifndef _circularProbe_hpp
#define _circularProbe_hpp

#include "probe.hpp"

//! \brief Circular Probe
//! \author V. Richefeu
class circularProbe : public Probe
{	
protected:
  
  double x_,y_,R_;
  
public:
  
  circularProbe(): x_(0.0), y_(0.0), R_(1.0) { }
  circularProbe(double x, double y, double R): x_(x), y_(y), R_(R) { }
  ~circularProbe() { }
   
  
  double area() {return (M_PI*R_*R_);}
  
  bool containDofCenter (dof * i)
  {
    double dx = i->mcx() - x_;
    double dy = i->mcy() - y_;
    return ( (dx*dx + dy*dy) < R_*R_);

  }
  
  bool containCenter (body2d * b)
  {
    double dx = b->x() - x_;
    double dy = b->y() - y_;
    return ( (dx*dx + dy*dy) < R_*R_);
  }
  
  bool contain (inter2d * k)
  {
    double dx = k->x() - x_;
    double dy = k->y() - y_;
    return ( (dx*dx + dy*dy) < R_*R_);
  }
  
  bool intersection ( body2d * b)
  {
	double dx = b->x() - x_;
    double dy = b->y() - y_;
	double d=dx*dx+dy*dy;
	//double r=b->sizeVerlet();
  
	return ( (R_-b->sizeVerlet())*(R_-b->sizeVerlet())< d  && d < (R_+b->sizeVerlet())*(R_+b->sizeVerlet())  );
	
//	return ( ( R_*R_< d || r*r< d) && d < R_*R_ + r*r + 2.*r*R_  );
  
  }
  
  bool containEntireBody ( body2d * b )
  {
	double dx = b->x() - x_;
    double dy = b->y() - y_;
  
	return ( dx*dx+dy*dy <  (b->sizeVerlet()-R_)*(b->sizeVerlet()-R_));
		//	&& dx*dx + dy*dy > pow( R_ + b->sizeVerlet(),2));
  }

bool probeInsideBody ( body2d * b )
  {
	double dx = b->x() - x_;
    double dy = b->y() - y_;

	return ( dx*dx+dy*dy  < (b->sizeVerlet()+R_)*(b->sizeVerlet()+R_) );
  }

  double & x()          { return x_; }
  double   x()    const { return x_; }
  double & y()          { return y_; }
  double   y()    const { return y_; }
  double & R()          { return R_; }
  double   R()    const { return R_; } 
	
};

#endif // _circularProbe_hpp
