#ifndef _dtk_point_hpp
#define _dtk_point_hpp

namespace dtk
{
  
  class point
  {
private :
    
    double x_;
    double y_;
    
public:
      
    point() : x_(0.0), y_(0.0) { }
	point(double x, double y) : x_(x),y_(y) {}
    
    double   x() const { return x_; }
    double & x()       { return x_; }
    double   y() const { return y_; }
    double & y()       { return y_; }       
  };
  
}

#endif // _dtk_point_hpp
