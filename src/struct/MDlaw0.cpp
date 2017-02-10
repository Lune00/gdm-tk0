#include "MDlaw0.hpp" 

void MDlaw0::checkRequiredParameters()
{
  
}

void MDlaw0::computeForces(unsigned int k)
{

}

/*
inline double MDlaw0::fn   (unsigned int k)
{
  unsigned int g1 = _nwk->inter(k)->first()->grp();
  unsigned int g2 = _nwk->inter(k)->second()->grp();
  double kn = _ref->getParameter("kn",g1,g2);
  double dn = _nwk->inter(k)->Dist();
  
  return -kn*dn;
}

inline double MDlaw0::ft   (unsigned int k)
{
  return 0.0;
}

inline double MDlaw0::frot (unsigned int k)
{
  return 0.0;
}	
*/
