#ifndef _heightProbe_hpp
#define _heightProbe_hpp

#include "probe.hpp"

//! \brief Probe defined by two horizontal lines of different height
//! \author C. Voivret
class heightProbe: public Probe
{	
protected:
  
  double h1_,h2_;
  //h1_ < h2_
  double width_;
    
public:
  
  heightProbe(): h1_(0.0), h2_(0.0),width_(0.0) { if (h2_<h1_) cout<<"@HeightProbe, Bad definiton of heightProbe"<<endl; }
  heightProbe(double x, double y): h1_(x), h2_(y){ if (h2_<h1_) cout<<"@HeightProbe, Bad definiton of heightProbe"<<endl; }
  heightProbe(double x, double y,double w): h1_(x), h2_(y), width_(w) { if (h2_<h1_) cout<<"@HeightProbe, Bad definiton of heightProbe"<<endl; }

  ~heightProbe() { }
  
  double & h1()          { return h1_; }
  double   h1()    const { return h1_; }
  double & h2()          { return h2_; }
  double   h2()    const { return h2_; }
  double & width()          { return width_; }
  double   width()    const { return width_; } 
  //double & area()          { return area_; }
  double   area()   { return (width_*(h2_-h1_)); } 
  
  double   halfHeight()    const { return .5*(h1_+h2_); }  
  
  bool containDofCenter (dof * i)
  {
    return ( (i->mcy() >= h1_ ) && ( i->mcy() <= h2_ ) );
  }
    
  bool containCenter (body2d * b)
  {
    
    return ( (b->y() >= h1_ ) && ( b->y() <= h2_ ) );
  }
  
  bool contain (inter2d * k)
  {
    
    return ( (k->y() >= h1_ ) && ( k->y() <= h2_ ) );
  }
  
  bool intersection ( body2d * b)
  {
	if( b->y()<h1_)
		return(  b->ymax() > h1_);
	else if ( b->y()>h2_)
		return(  b->ymin() < h2_);
	else 
		return( false );
  
  }
  
  bool containEntireBody ( body2d * b )
  {
	return ( ( b->ymin() <= h1_) && (b->ymax() <= h2_) );
	}


};

#endif // _heightProbe_hpp

