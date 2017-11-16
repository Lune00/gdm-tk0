#ifndef _inter2d_h
#define _inter2d_h

#include <string>
#include "body2d.hpp"
#include "dof.hpp"
#include <typeinfo>
//#include "groupRelationData.hpp"
//#include "dof.hpp"
// definition d'identificateur de type int remplacement du type mentionné par un int
#define _type_dkdk 0
#define _type_dkdkP 1
#define _type_dkrl 2
#define _type_pgpg 3
#define _type_rlrl 4
#define _type_pgrl 5
#define _type_dkdkco 6

class GroupRelationData;

class inter2d
{

	protected:

		unsigned int rank_;    //!< Number of active contacts
		unsigned int current_; //!< Current contact point (0 or 1)
		bool activated_; //!< An activated inter have 0 (distance interaction), 1 or 2 contacts
		double dAct_; //!< Distance of activation (0 by default)
		double longeur_; // Longeur de la surface de contact

		unsigned int type_;
		double x_,y_;
		double nx_,ny_;
		double tx_,ty_;
		double fn_,ft_,frot_;
		double vn_,vt_,vrot_;

	public:

		virtual double  Dist()     = 0;
		virtual double  Overlap()  = 0; // return positive value if overlap  else -1
		virtual void    Frame()    = 0;
		virtual void    Frame2()    = 0;
		virtual bool    Activate()    = 0;
		virtual void    Kin()      = 0;
		virtual void    Res()      = 0;
		virtual void    Res(const double, const double, const double = 0.0) = 0;
		virtual void    CDcoeff()  = 0;
		virtual void    CDcoeff(GroupRelationData*) {} ;
		virtual double  An(const double) = 0;
		virtual double  At(const double) = 0;
		virtual double  As(const double) = 0;

		virtual body2d* first()      = 0;
		virtual body2d* second()     = 0;

		virtual void plug(body2d* b1, body2d* b2) = 0;

		virtual body2d* lexifirst()  = 0;
		virtual body2d* lexisecond() = 0;

		virtual void read(istream & is, unsigned int * Id1, unsigned int * Id2) = 0;
		virtual void write(ostream & os)    = 0;
		virtual void writeMGP(ostream & os) = 0;

		virtual void clearForceAndMoment() = 0; 

		virtual double & x()        = 0;
		virtual double   x()  const = 0;
		virtual double & y()        = 0;
		virtual double   y()  const = 0;

		virtual double & vn()         = 0;
		virtual double   vn()   const = 0;
		virtual double & vt()         = 0;
		virtual double   vt()   const = 0;
		double & vrot()       { return vrot_; }
		double   vrot() const { return vrot_; }
		virtual double & fn()         = 0;
		virtual double   fn()   const = 0;
		virtual double & ft()         = 0;
		virtual double   ft()   const = 0;
		double & frot()       { return frot_; }
		double   frot() const { return frot_; }
		virtual double   Vbranchx()   const = 0;
		virtual double   Vbranchy()   const = 0;

		virtual double   fx()   const= 0;
		virtual double   fy()   const= 0; 

		virtual ~inter2d() { }  

		inter2d()
		{
			fn_ = ft_ = frot_ = 0.0;
			vn_ = vt_ = vrot_ = 0.0;
			rank_ = 0;
			longeur_ = 0;
			activated_=false;
			dAct_=0.;
			type_ = 9;
		}

		static inter2d* factory(body2d* b1, body2d* b2);
		static inter2d* factoryP(body2d* b1, body2d* b2,double P);
		static inter2d* factory(string type);

		unsigned int   rang() const { return rank_; }
		unsigned int & rang()       { return rank_; }

		unsigned int   current() const { return current_; }
		unsigned int & current()       { return current_; }

		double   longeur() const {return longeur_;}
		double & longeur()       {return longeur_;}

		bool   activate() const { return activated_; }
		bool & activate()       { return activated_; }

		double & dAct()       { return dAct_; }
		double   dAct() const { return dAct_; }

		double & nx()       { return nx_; }
		double   nx() const { return nx_; }
		double & ny()       { return ny_; }
		double   ny() const { return ny_; }
		double & tx()       { return tx_; }
		double   tx() const { return tx_; }
		double & ty()       { return ty_; }
		double   ty() const { return ty_; }
		unsigned int & type()       { return type_;}
		unsigned int   type() const { return type_;}

		struct compareLength 
		{ 
			bool operator () (const inter2d & a1, const inter2d & a2) const 
			{ 
				return (sqrt(a1.Vbranchx()*a1.Vbranchx()+ a1.Vbranchy()*a1.Vbranchy()) < 
						sqrt(a2.Vbranchx()*a2.Vbranchx()+ a2.Vbranchy()*a2.Vbranchy()));
			}

			bool operator () (const inter2d * a1, const inter2d * a2) const 
			{ 
				return (sqrt(a1->Vbranchx()*a1->Vbranchx()+a1->Vbranchy()*a1->Vbranchy()) < 
						sqrt(a2->Vbranchx()*a2->Vbranchx()+a2->Vbranchy()*a2->Vbranchy()));
			}	
		};

		struct compareNormalForce 
		{ 
			bool operator () (const inter2d & a1, const inter2d & a2) const 
			{ 
				return (a1.fn() < a2.fn());
			}

			bool operator () (const inter2d * a1, const inter2d * a2) const 
			{ 
				return (a1->fn() < a2->fn());
			}	
		};



};

#endif // _inter2d_h

