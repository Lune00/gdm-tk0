#ifndef _control_hpp
#define _control_hpp

#include <iostream>
#include <string>
#include "talk.hpp"

#define _FORCE      0
#define _VELOCITY   1
#define _PRESSURE   2

using namespace std;

//! \brief Define the degrees of freedom for one body
//! \author V. Richefeu
class control
{
  
public:
  control(unsigned int x,unsigned int y, double xval, double yval) :
  x_(x), y_(y), xval_(xval), yval_(yval) { rot_ = _VELOCITY; rotval_ = 0.0; }
  
  control(unsigned int x,unsigned int y,unsigned int rot,
          double xval, double yval, double rotval) :
  x_(x), y_(y), rot_(rot), xval_(xval), yval_(yval), rotval_(rotval) { }
  
  control() 
  { 
    x_    = y_    = rot_    = _VELOCITY;
    xval_ = yval_ = rotval_ = 0.0; 
  }
  
  unsigned int      x() const { return x_;    }
  unsigned int    & x()       { return x_;    }
  unsigned int      y() const { return y_;    }
  unsigned int    & y()       { return y_;    }
  unsigned int    rot() const { return rot_;  }
  unsigned int  & rot()       { return rot_;  }

  double         xval() const { return xval_; }
  double       & xval()       { return xval_; }
  double         yval() const { return yval_; }
  double       & yval()       { return yval_; }
  double       rotval() const { return rotval_; }
  double     & rotval()       { return rotval_; }
  
  // Some usefull funtions for debugging
  void check();
  void print();
    
private:
    
  unsigned int x_    , y_    , rot_;
  double       xval_ , yval_ , rotval_;
};

#endif // _control_hpp
