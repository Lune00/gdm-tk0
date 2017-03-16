#ifndef _interdof_h
#define _interdof_h
#include "dof.hpp"
#include "inter2d.hpp"
#include "body2d.hpp"
#include "system.hpp"

class interdof
{
	dof * i_;
	dof * j_;
   
	vector<inter2d * > linter_;
	
	double nx_,ny_, tx_, ty_;
	double fresNormal_;
	double fresOrthoNorm_;
	double fresNormalx_;
	double fresNormaly_;
	double vbranchx_;
	double vbranchy_;
	double vbranch_;
  
  public:
  
  interdof(dof * i,dof * j):i_(i), j_(j) {fresNormal_ = fresOrthoNorm_ = vbranchx_ = vbranchy_ = fresNormalx_ = fresNormaly_ = 0.;}
  dof* first()  { return i_; } // const
  dof* second() { return j_; }
  
  vector <inter2d *> & linter() {return linter_; }
  int rank() {return linter_.size();}

  void calcFres();
  void Frame();
  double fresNormal(){ return fresNormal_;}
  double fresOrthoNorm(){return fresOrthoNorm_;}
  double nx() const {return nx_;}
  double ny() const {return ny_;}
  double tx() const {return -ny_;}
  double ty() const {return nx_;}
  double vbranchx()  {return vbranchx_;}
  double vbranchy()  {return vbranchy_;}
  double fresNormalx() {return fresNormalx_;}
  double fresNormaly() {return fresNormaly_;}
  double vbranch()	 {return vbranch_;}
};

#endif 


