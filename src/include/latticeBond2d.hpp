#ifndef _latticeBond2d_h
#define _latticeBond2d_h

#include <iostream>
#include <string>

using namespace std;

class latticeBond2d
{

 protected:

  unsigned int _id;
  unsigned int _type;
  unsigned int _upScaleId;
  latticeNode2d *_first,*_second;
  double _K,_F;
  double _F0;

 public:

  void read(istream & is)     = 0;
  void write(ostream & os)    = 0;
  void writeMGP(ostream & os) = 0; 
  void writePS(ostream & os)  = 0;
	
  ~latticeBond2d() { }
   latticeBond2d() { }

  unsigned int & id()              { return _id; }	
  unsigned int   id()        const { return _id; }
  unsigned int & type()            { return _type; }	
  unsigned int   type()      const { return _type; }
  unsigned int & upScaleId()       { return _upScaleId; }	
  unsigned int   upScaleId() const { return _upScaleId; }  
	
  double & x()          { return _x; }
  double   x()    const { return _x; }
  double & y()          { return _y; }
  double   y()    const { return _y; }

};

#endif // _latticeBond2d_h



