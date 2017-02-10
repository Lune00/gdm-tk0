#ifndef _latticealgo_h
#define _latticealgo_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "lattice.hpp"

class LatticeAlgo
{	
  Lattice*  _lat;

  double        _dt,_tm;
  unsigned int  _nitermn,_nitermx,_niter;
  unsigned int  _nsi,_nsf,_ns;
		
public:

  Latticealgo(Lattice* lat, double dt) : _lat(lat), _dt(dt)
  { _tm=0.0; }
	
  void read_parameters(istream&);

  void         stand();
  void         look();
  void         hand();
  void         step(unsigned int);
	
  unsigned int contactStatut( inter2d *);

  unsigned int minimizeEp();


  unsigned int    nsi()     const { return _nsi; }
  unsigned int  & nsi()           { return _nsi; }
  unsigned int    ns()      const { return _ns; }
  unsigned int  & ns()            { return _ns; }
  unsigned int    nsf()     const { return _nsf; }
  unsigned int  & nsf()           { return _nsf; }
  double          dt()      const { return _dt; }
  double        & dt()            { return _dt; }
  unsigned int    nitermn() const { return _nitermn; }
  unsigned int  & nitermn()       { return _nitermn; }
  unsigned int    nitermx() const { return _nitermx; }
  unsigned int  & nitermx()       { return _nitermx; }
  unsigned int    niter()   const { return _niter; }
  unsigned int  & niter()         { return _niter; }
  double          tm()      const { return _tm; }
  double        & tm()            { return _tm; }
	
};

#endif // _CDalgo_h
