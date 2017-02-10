#ifndef _dof_hpp
#define _dof_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include "talk.hpp"
#include "body2d.hpp"
#include "inter2d.hpp"
#include "vertex.hpp"

#define _FORCE      0
#define _VELOCITY   1
#define _PRESSURE   2


using namespace std;

//! \brief Define the degrees of freedom for one or a group of bodies
//! \author C. Voivret
//class body2d;
class dof
{

private:
//Structure similaire au "control" pour imposer une force ou un deplacement
	unsigned int x_    , y_    , rot_;
	double       xval_ , yval_ , rotval_;

	double m_,mom_,area_;
	double xGravity_,yGravity_;
	unsigned int id_;

//Coordonnees du centre de masse
	double mcx_,mcy_,mcrot_;//x y rot

//vitesses communes au corps controlés
	double mcvx_,mcvy_,mcvrot_;
//Forces appilquées au centre de masse du dof	
	double mcfx_,mcfy_,mcfrot_;

	vector <body2d*> ctrlBodies_;
	
	vector<gdm::vertex> vertex_;
//dof periodique ? :
// - oui : rotation automatiquement bloquée, le corps va d'un bout a l'autre de la periode
// - non : les rotations sont permises et le corps est periodiquement translater dans son ensemble
	bool isPeriodic_;
	
	bool multiBody_;


public:
	dof(unsigned int x,unsigned int y, double xval, double yval) :
	x_(x), y_(y), xval_(xval), yval_(yval) 
		{ rot_ = _VELOCITY; rotval_ = 0.0;m_=1.;mcfx_=0.;mcfy_=0;mcfrot_=0;xGravity_= 0.;yGravity_=0.;multiBody_=true;area_=-1; }

	dof(unsigned int x,unsigned int y,unsigned int rot,
		double xval, double yval, double rotval) :
	x_(x), y_(y), rot_(rot), xval_(xval), yval_(yval), rotval_(rotval) {m_=1.;mcfx_=0.;mcfy_=0;mcfrot_=0;xGravity_= 0.;yGravity_=0.;multiBody_=true;area_=-1; }

	dof() 
	{ 
		x_    = y_    = rot_    = _FORCE;
		xval_ = yval_ = rotval_ = 0.0;
		isPeriodic_=false;
		multiBody_=true;
		area_=-1;
	}
//Structure de controle
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
// Forces vitesses et position du centre de masse
	double         mcfx()   const { return mcfx_; }
	double       & mcfx()         { return mcfx_; }
	double         mcfy()   const { return mcfy_; }
	double       & mcfy()         { return mcfy_; }
	double         mcfrot() const { return mcfrot_; }
	double       & mcfrot()       { return mcfrot_; }
	
	double         mcvx()   const { return mcvx_; }
	double       & mcvx()         { return mcvx_; }
	double         mcvy()   const { return mcvy_; }
	double       & mcvy()         { return mcvy_; }
	double         mcvrot() const { return mcvrot_; }
	double       & mcvrot()       { return mcvrot_; }
	
	double      mcx() const { return mcx_;    }
	double    & mcx()       { return mcx_;    }
	double      mcy() const { return mcy_;    }
	double    & mcy()       { return mcy_;    }
	double      mcrot() const { return mcrot_;    }
	double    & mcrot()       { return mcrot_;    }
//masse et moment d'inertie
	double      m() const { return m_;    }
	double    & m()       { return m_;    }
	double      mom() const { return mom_;    }
	double    & mom()       { return mom_;    }
	double      Area() const { return area_;    }
	double    & Area()       { return area_;    }
	

	unsigned int    id() const { return id_;  }
	unsigned int  & id()       { return id_;  }

	double & xG()       { return xGravity_;}
	double   xG() const { return xGravity_;}

	double & yG()       { return yGravity_;}
	double   yG() const { return yGravity_;}
	
	bool & isPeriodic()       { return isPeriodic_;}
	bool   isPeriodic() const { return isPeriodic_;}
	
	bool & multiBody()       { return multiBody_;}
	bool   multiBody() const { return multiBody_;}
	

	vector<body2d*> & lctrlBodies()       { return ctrlBodies_; } 
	vector<body2d*>   lctrlBodies() const { return ctrlBodies_; }
	body2d* ctrlBody(unsigned int i) { return ctrlBodies_[i]; }

	body2d* lowerBody();
	
	void affect(unsigned int x,unsigned int y,unsigned int rot,
		double xval, double yval, double rotval);

	void resDof(double, double, double, double ,double ); //x et y interaction
	void move( double );
	void imposeForce();
	void imposeVelocity();
	void imposeVelocityOfForce(double );
	void imposeForceOfVelocity();

    void ComputeImposedMass(double);
    
	void plugBody( unsigned int ,string ,vector <unsigned int> &);//A faire avec des flux....
	void exportBodyId( string );
	void plugBody(body2d *);

	void setGravity( double xg,double yg) { xGravity_=xg;yGravity_=yg;}
	void fill( double density);
	void computeMassCenter();
	void computeMoment();
	void computeVertex();

// Corps les plus a gauche et a droite 
	body2d* rightBody();
	body2d* leftBody();
	void translate( double);
	void translateBodies( double, double,double );
	
// Some usefull functions for debugging
	void check();
	void print();
	
//  Ecriture ds spl
	unsigned int writeSpl( ostream & );
};

#endif // _dof_hpp

