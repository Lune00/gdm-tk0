#ifndef _rline_h
#define _rline_h

#include <math.h> 
#include "body2d.hpp"

class rline : public body2d
{

  double R_;
  double L_; // Longueur tot ou demi longueur ??? a verifier...
		
  rline() : body2d(), R_(1.0), L_(0.0) {type()=_type_rline; }
  friend class body2d;

public:
	rline( double r, double l) : body2d(), R_(r), L_(l) { type()=_type_rline;}

  ~rline() {}
  //rline(body2d &) { }

  double   R() const { return R_; }
  double & R()       { return R_; }
  double   L() const { return L_; }
  double & L()       { return L_; }  
	
  void read    (istream & is);
  void write   (ostream & os);
  void writeMGP(ostream & os);
  void writePS (ostream & os);
    void writeM (ostream & os);
  rline* duplicate();
  
  double Area() const { return M_PI * R_ * R_ + 2.0 * L_ * R_; } 
 
  double xmin() { return x_ - L_ * fabs(cos(rot_)) - R_; }
  double xmax() { return x_ + L_ * fabs(cos(rot_)) + R_; }
  double ymin() { return y_ - L_ * fabs(sin(rot_)) - R_; }
  double ymax() { return y_ + L_ * fabs(sin(rot_)) + R_; }
  double sizeVerlet() const { return R_; }
	
  void Fill(double);	
};

#endif // _rline_h
