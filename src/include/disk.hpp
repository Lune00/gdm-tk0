#ifndef _disk_hpp
#define _disk_hpp

#include <math.h>
#include "body2d.hpp"

//! \brief A disk shaped body
//! \author V. Richefeu
class disk : public body2d
{
  
private:
  
  double R_; //< Radius of the disk
 
  disk(body2d &) {type()=_type_disk; }
  disk() : body2d(), R_(1.0) {type()=_type_disk; }
  friend class body2d;

public:

  ~disk() {}
	disk( double x_ , double y_, double r_): body2d(),R_(r_){ x()=x_;y()=y_;}

 double   R() const { return R_; }
/*
 
  double   x() const {return x_;}
  double   y() const {return y_;}
  double   rot() const {return rot_;}
  double   vrot() const {return vrot_;}
  double   vx() const {return vx_;}
  double   vy() const {return vy_;}
  */
  
  double & R()       { return R_; }
	
  void read(istream & is);
  void write(ostream & os);
  void writeMGP(ostream & os);
  void writePS(ostream & os); 
  void writeM(ostream & os);
  disk* duplicate();
  double sizeVerlet() const { return R_; }
 
  double Area() const { return M_PI * R_ * R_; } 
  
  double xmin()  { return x_ - R_; }
  double xmax()  { return x_ + R_; }
  double ymin()  { return y_ - R_; }
  double ymax()  { return y_ + R_; }

 



  void Fill(double);

struct compareRadius 
    { 
      bool operator () (const disk & a1, const disk & a2) const 
      { 
        return (a1.R() < a2.R());
      }

      bool operator () (const disk * a1, const disk * a2) const 
      { 
        return (a1->R() < a2->R());
      }	
    
};
};

#endif // _disk_hpp
