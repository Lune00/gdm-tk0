#ifndef _latticeNode2d_h
#define _latticeNode2d_h

#include <iostream>
#include <string>

using namespace std;

class latticeNode2d
{

 protected:

  unsigned int _id;
  unsigned int _type;
  unsigned int _upScaleId;
  double _x,_y;
  double _Ep;

	
 public:

  void read(istream & is)     = 0;
  void write(ostream & os)    = 0;
  void writeMGP(ostream & os) = 0; 
  void writePS(ostream & os)  = 0;
	
  ~latticeNode2d()  { }
   latticeNode2d()  { }

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

#endif // _latticeNode2d_h



