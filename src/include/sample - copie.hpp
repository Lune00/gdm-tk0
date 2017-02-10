#ifndef _sample_hpp
#define _sample_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include "body2d.hpp"
#include "dof.hpp"

#include "disk.hpp" // pourquoi ?? le display peut etre (a mettre au clair)
#include "polyg.hpp" // ?
#include "rline.hpp" // ?
#include "purge.h"

//! \brief A object containing the list of bodies
//! \author V. Richefeu
class Sample
{
  vector<body2d*> lbody_;
  vector<unsigned int> leftband_;
  vector<unsigned int> rightband_;
  double xmin_, xmax_, ymin_, ymax_;
  double rmin_,rmax_,rmoy_;

  
  bool   isMonoPeriodic_;
  double leftBoundary_;
  double rightBoundary_;
  double boundWidth_;    //< Largeur de la periode
  double bandWidth_;     //< Largeur de la bande de detection
  
public:
    
  Sample(): isMonoPeriodic_(false) { }
  ~Sample() { purge(lbody_); }
  
  vector<body2d*> & lbody()       { return lbody_; } 
  vector<body2d*>   lbody() const { return lbody_; }
  body2d* body(unsigned int i) { return lbody_[i]; }
  
  void substituteBody(unsigned int i, body2d& B) 
    {
    if (i >= lbody_.size()) return;
    delete lbody_[i];
    lbody_[i] = &B;
    }
  
  vector<unsigned int> & leftband()       { return leftband_; } 
  vector<unsigned int>   leftband() const { return leftband_; }
  unsigned int leftband( unsigned int i)  { return leftband_[i];}
  
  vector<unsigned int> & rightband()       { return rightband_; } 
  vector<unsigned int>   rightband() const { return rightband_; }
  unsigned int rightband( unsigned int i)  { return rightband_[i];}
	
  vector <dof*> read  (istream&);
  void write (ostream&);
  void write (const char * fname);
  void fill  (double);
  
  //! \brief Compute the boundaried (xmin, xmax, ymin, ymax) of the sample
  void updateBoundaries();
  
  //! \brief Compute Rmin and Rmax
  //! \author Charles Voivret 
  void radiusExtrema(unsigned int );
  
  void updateBands();
  void translate();
  
  void swapBodies(unsigned int,unsigned int);
  
  Sample(double left,double right,double bandW): isMonoPeriodic_(true),leftBoundary_(left),
    rightBoundary_(right),boundWidth_(right-left),bandWidth_(bandW){ }
  
  //for sample construct with geometrical rules
  //! \author Charles Voivret
  void definePeriodicityCV(bool periodic);
  
  void SortByHeight();
  
  double   leftBoundary()   const { return leftBoundary_; }
  double & leftBoundary()         { return leftBoundary_; }
  
  double   rightBoundary()  const { return rightBoundary_; }
  double & rightBoundary()        { return rightBoundary_; }
  
  double   bandWidth()      const { return bandWidth_; }
  double & bandWidth()            { return bandWidth_; }
  
  double   boundWidth()     const { return boundWidth_; }
  double & boundWidth()           { return boundWidth_; }
  
  double   rmax()           const { return rmax_;}
  double & rmax()                 { return rmax_;}
  
  double   rmin()           const { return rmin_;}
  double & rmin()                 { return rmin_;}

double   rmoy()           const { return rmoy_;}
  double & rmoy()                 { return rmoy_;}

  
  bool   isMonoPeriodic()   const { return isMonoPeriodic_; }
  bool & isMonoPeriodic()         { return isMonoPeriodic_; }
  
	void addBody( body2d * bod)
	 { lbody_.push_back(bod);
		//cout<<"taille lbody "<<lbody_.size()<<endl;
       //lbody_[n].back()->id() = ;
		}
		void removeBody( body2d * );
  
  double  xmin() const { return xmin_; }
  double  xmax() const { return xmax_; }
  double  ymin() const { return ymin_; }
  double  ymax() const { return ymax_; }	
};


#endif // _sample_hpp
