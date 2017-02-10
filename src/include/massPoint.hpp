#ifndef _massPoint_hpp
#define _massPoint_hpp

#include <math.h>
#include "body2d.hpp"

//! \brief A point with a mass
//! \author V. Richefeu
class massPoint : public body2d
{
  
private:
 
  double R_; //< an imaginary radius
  
  massPoint(body2d &) { }
  massPoint() : body2d(), R_(1.0) { }

  friend class body2d;

public:

  ~massPoint() {}

  double   R() const { return R_; }
  double & R()       { return R_; }
 double sizeVerlet() const { return R_; }

  void read(istream & is);
  void write(ostream & os);
  void writeMGP(ostream & os);
  void writePS(ostream & os); 
    void writeM (ostream & os);
  massPoint* duplicate();
  
  double Area() const { return M_PI * R_ * R_; } 
  /*double xmin() const { return x_ - R_; }
  double xmax() const { return x_ + R_; }
  double ymin() const { return y_ - R_; }
  double ymax() const { return y_ + R_; }*/

  double xmin() { return x_ - R_; }
  double xmax() { return x_ + R_; }
  double ymin() { return y_ - R_; }
  double ymax() { return y_ + R_; }
  void Fill(double);
};

#endif // _massPoint_hpp
