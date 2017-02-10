#ifndef _gdm_vertex_hpp
#define _gdm_vertex_hpp

#include <math.h>

namespace gdm
{
  
  class vertex
  {
private :
    
    double _x;
    double _y;
    
public:
      
    vertex() : _x(0.0), _y(0.0) { }
    vertex(double x, double y) : _x(x),_y(y) {}
    void Normalize()
      {
      }
	double Norm()
	{
		return sqrt( _x*_x + _y*_y);
	}
    
    double   x() const { return _x; }
    double & x()       { return _x; }
    double   y() const { return _y; }
    double & y()       { return _y; }    
      
  };
  
}

#endif // _gdm_vertex_hpp


