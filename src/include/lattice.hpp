#ifndef _lattice_h
#define _lattice_h

#include <iostream>
#include <fstream>
#include <vector>
#include "latticeNode2d.hpp"
#include "latticeBond2d.hpp"
#include "purge.h"

class Lattice
{
  vector<latticeNode2d*> _lnode;
  vector<latticeBond2d*> _lbond;
  double _xmin,_xmax,_ymin,_ymax;
	
 public:

  Lattice() { }
  ~Lattice() { purge(_lnode); purge(_lbond); }

  vector<latticeNode2d*> & lnode()       { return _lnode; } 
  vector<latticeNode2d*>   lnode() const { return _lnode; }
  latticeNode2d* node(unsigned int i) { return _lnode[i]; }
  
  vector<latticeBond2d*> & lbond()       { return _lbond; } 
  vector<latticeBond2d*>   lbond() const { return _lbond; }
  latticeBond2d* bond(unsigned int i) { return _lbond[i]; }  
	
  void read(istream&);
  void write(const char * fname);
  void init();
  void updateBoundaries();

  double  xmin() const { return _xmin; }
  double  xmax() const { return _xmax; }
  double  ymin() const { return _ymin; }
  double  ymax() const { return _ymax; }	
};


#endif // _lattice_h
