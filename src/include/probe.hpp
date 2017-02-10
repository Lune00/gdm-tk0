#ifndef _probe_hpp
#define _probe_hpp

#include "body2d.hpp"
#include "inter2d.hpp"
#include "dof.hpp"

//! \brief Virtual class for probes
//! \author V. Richefeu
class Probe
{	
  
public:
  
  Probe() { }
  virtual ~Probe() { }

  virtual bool containCenter (body2d *) = 0;	
  virtual bool contain (inter2d *) = 0;
  virtual bool containDofCenter (dof *) =0;	
  virtual double area ()=0;
  
};

#endif // _probe_hpp
