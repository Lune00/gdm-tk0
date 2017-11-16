#ifndef _polyg_hpp
#define _polyg_hpp

//#include <math.h>
#include <cmath>
#include <vector>
#include "body2d.hpp"
#include "vertex.hpp"

//! \brief A polygonal shaped body
//! \author V. Richefeu
class polyg : public body2d
{

private:
    
  vector<gdm::vertex> vertex_; //!< Body shape is this vertexes ({0,0} is the mass center) 
  vector<gdm::vertex> normal_; //!< Outgoing normals ; vecto don vi huong ra ngoai thu i
  double Rout_;                //!< Circumscribed radius
  double Rin_;                //!< Circumscribed radius
	
  double Rmass_;               //!< Radius defined for a disk equivalent mass
  double Area_;
  double xmin_,xmax_,ymin_,ymax_; //Ajuster le 13/1/2012 par NGUYEN
  double c,s;

  
  friend class body2d;
	
public:
  polyg(body2d &) { type()=_type_polyg;}
  polyg() : body2d() {type()=_type_polyg; }


  ~polyg() {}
  
  double   Rout() const { return Rout_; }
  double & Rout()       { return Rout_; }

  vector<gdm::vertex> & Vertex()       { return vertex_; }
  vector<gdm::vertex>   Vertex() const { return vertex_; }
  
  gdm::vertex & Vertex(unsigned int i)       { return vertex_[i]; }
  gdm::vertex   Vertex(unsigned int i) const { return vertex_[i]; }

  vector<gdm::vertex> & Normal()       { return normal_;    }
  vector<gdm::vertex>   Normal() const { return normal_;    }
  
  gdm::vertex & Normal(unsigned int i)       { return normal_[i]; }
  gdm::vertex   Normal(unsigned int i) const { return normal_[i]; }
  
  void read(istream & is);
  void write(ostream & os);
  void writeMGP(ostream & os);
  void writePS(ostream & os);
    void writeM(ostream & os);
  
  polyg* duplicate();
  
  double Area() const { return M_PI * Rmass_ * Rmass_; } 
  /*double xmin() const { return x_ - Rout_; }
  double xmax() const { return x_ + Rout_; }
  double ymin() const { return y_ - Rout_; }
  double ymax() const { return y_ + Rout_; }*/



  double xmin() 
  { 
  	c=cos(rot_),s=sin(rot_);
  	xmin_=x_+vertex_[0].x()*c-vertex_[0].y()*s;
  	for (unsigned int i=1;i<vertex_.size();i++)
  	{
  		xmin_=(xmin_>x_+vertex_[i].x()*c-vertex_[i].y()*s) ? x_+vertex_[i].x()*c-vertex_[i].y()*s : xmin_;
  	}
  	return xmin_;
  }
  
  double xmax()  
  {
  	c=cos(rot_),s=sin(rot_);
    xmax_=x_+vertex_[0].x()*c-vertex_[0].y()*s;
  	for (unsigned int i=1;i<vertex_.size();i++)
  	{
  		xmax_=(xmax_<x_+vertex_[i].x()*c-vertex_[i].y()*s) ? x_+vertex_[i].x()*c-vertex_[i].y()*s : xmax_;
  	}
  	return xmax_;
  }
  
  double ymin() 
  { 
    c=cos(rot_),s=sin(rot_);
  	ymin_=y_+vertex_[0].x()*s+vertex_[0].y()*c;
  	for (unsigned int i=1;i<vertex_.size();i++)
  	{
  		ymin_=(ymin_>y_+vertex_[i].x()*s+vertex_[i].y()*c) ? y_+vertex_[i].x()*s+vertex_[i].y()*c : ymin_;
  	}
  	return ymin_;
  } 
  
  double ymax() 
  { 
    c=cos(rot_),s=sin(rot_);
  	ymax_=y_+vertex_[0].x()*s+vertex_[0].y()*c;
  	for (unsigned int i=1;i<vertex_.size();i++)
  	{
  		ymax_=(ymax_<y_+vertex_[i].x()*s+vertex_[i].y()*c) ? y_+vertex_[i].x()*s+vertex_[i].y()*c : ymax_;
  	}
  	return ymax_;
  }
  
  double sizeVerlet() const { return Rout_; }
	
  void Fill(double);
  
  void adjustCenter ();
  void convexify();
  
  bool Contain(gdm::vertex &v);
  vector<gdm::vertex> SectionSegment(gdm::vertex d1,gdm::vertex d2);
  void Arrangement();
  int ispolysimple(); 
  void Invert();
};

#endif // _polyg_hpp
