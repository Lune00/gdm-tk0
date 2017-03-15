#ifndef _body2d_hpp
#define _body2d_hpp

#include <iostream>
#include <string>
using namespace std;

#define _type_disk 0
#define _type_polyg 1
#define _type_rline 2

//! \brief Virtual class for two-dimensionnal bodies
//! \author V. Richefeu
class dof;
class body2d
{

protected:

	unsigned int grp_;     //< Group identifier
	unsigned int id_;      //< Body identifier
	double x0_,y0_,rot0_;
	double x_,y_,rot_;
	double vx_,vy_,vrot_;
	double ax_,ay_,arot_;
	double fx_,fy_,frot_;
	double mass_,mom_;
	unsigned int type_;
	unsigned int z_; // nombre de contact, seulement pour analyse (non pour algo)
	dof* dof_;

public:

	virtual double Area() const = 0;
	/*virtual double xmin() const = 0;
	virtual double xmax() const = 0;
	virtual double ymin() const = 0;
	virtual double ymax() const = 0;*/

	virtual double xmin() = 0;
	virtual double xmax() = 0;
	virtual double ymin() = 0;
	virtual double ymax() = 0;
//! \brief Return the outer size 
  	virtual double sizeVerlet() const = 0;
	
	//virtual  vector<gdm::vertex> & Vertex()       { return vertex_; }
	//virtual  gdm::vertex & Vertex(unsigned int i)       { return vertex_[i]; }
	
//! \brief Compute the mass and moment of the body
	virtual void   Fill(double) = 0;

	virtual void read(istream & is)     = 0;
	virtual void write(ostream & os)    = 0;
	virtual void writeMGP(ostream & os) = 0; 
	virtual void writePS(ostream & os)  = 0;
    virtual void writeM(ostream & os) = 0;
	virtual body2d* duplicate() = 0;

	virtual ~body2d() {}

	body2d()
	{
		x0_ = y0_ = rot0_ = 0.0; 
		vx_ = vy_ = vrot_ = 0.0;
		fx_ = fy_ = frot_ = 0.0;
		mass_ = mom_ = 1.0;
		dof_=0;
		z_ = 0 ;
	}

	static body2d* factory(string type);

	struct compareHeight 
	{ 
		bool operator () (const body2d & a1, const body2d & a2) const 
		{ 
			return (a1.y() < a2.y());
		}

		bool operator () (const body2d * a1, const body2d * a2) const 
		{ 
			return (a1->y() < a2->y());
		}	
	}; 

	unsigned int & grp()       { return grp_; }	
	unsigned int   grp() const { return grp_; }

	unsigned int & id()       { return id_; }	
	unsigned int   id() const { return id_; }

	double & mass()       { return mass_; }
	double   mass() const { return mass_; }	
	double & mom()        { return mom_;  }
	double   mom()  const { return mom_;  }

	double & x0()          { return x0_; }
	double   x0()    const { return x0_; }
	double & y0()          { return y0_; }
	double   y0()    const { return y0_; }
	double & rot0()        { return rot0_; }
	double   rot0()  const { return rot0_; }

	double & x()          { return x_; }
	double   x()    const { return x_; }
	double & y()          { return y_; }
	double   y()    const { return y_; }
	double & rot()        { return rot_; }
	double   rot()  const { return rot_; }
	double & fx()         { return fx_; }
	double   fx()   const { return fx_; }
	double & fy()         { return fy_; }
	double   fy()   const { return fy_; }
	double & frot()       { return frot_; }
	double   frot() const { return frot_; }
	double & ax()         { return ax_; }
	double   ax()   const { return ax_; }
	double & ay()         { return ay_; }
	double   ay()   const { return ay_; }
	double & arot()       { return arot_; }
	double   arot() const { return arot_; }
	double & vx()         { return vx_; }
	double   vx()   const { return vx_; }
	double & vy()         { return vy_; }
	double   vy()   const { return vy_; }
	double & vrot()       { return vrot_; }
	double   vrot() const { return vrot_; }
	dof *  & bodyDof()       { return dof_; }
	dof *    bodyDof() const { return dof_; }
    unsigned int & type()    { return type_;}
    unsigned int   type() const { return type_;}

    unsigned int & z() {return z_;}
    unsigned int z() const {return z_;}
};

#endif // _body2d_hpp



