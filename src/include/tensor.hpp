#ifndef _tensor_hpp
#define _tensor_hpp
#include <math.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace gdm
{
  
  //! \brief A 2x2 Tensor
  //! \author V. Richefeu
  class Tensor2x2
  {
private :
    
    string name_;
    double xx_, xy_;
    double yx_, yy_;
	double l1_, l2_;
	double vx1_, vy1_;
	double vx2_, vy2_;
    
public:
      
	Tensor2x2() : xx_(0.0), xy_(0.0), yx_(0.0), yy_(0.0), l1_(0.0), l2_(0.0){ name_ = "noName";}
    Tensor2x2(string name) : name_(name), xx_(0.0), xy_(0.0), yx_(0.0), yy_(0.0), l1_(0.0), l2_(0.0) { }
//	~Tensor2x2()
//		{cout<<" destruction de " <<name_<<endl;}
//	void adresse() {cout<<"adresse de "<<name_<<" "<<this<<endl;};
    void print();
	void symmetrize(); 
	void scalarMult( double );   
	void eigenValues( );
	void eigenVectors( );
	double  majorDirection();
	double  Tr(){ return (xx_+yy_);}
    
    double   xx() const { return xx_; }
    double & xx()       { return xx_; }
    
    double   xy() const { return xy_; }
    double & xy()       { return xy_; }
    
    double   yx() const { return yx_; }
    double & yx()       { return yx_; }
    
    double   yy() const { return yy_; }
    double & yy()       { return yy_; }
	
	double   l1() const { return l1_; }
    double & l1()       { return l1_; }
	
	double   l2() const { return l2_; }
    double & l2()       { return l2_; }
	
	double   vx1() const { return vx1_; }
    double & vx1()       { return vx1_; }
	
	double   vx2() const { return vx2_; }
    double & vx2()       { return vx2_; }
	
	double   vy1() const { return vy1_; }
    double & vy1()       { return vy1_; }
	
	double   vy2() const { return vy2_; }
    double & vy2()       { return vy2_; }

	
	string   name() const { return name_; }
    string & name()       { return name_; }
	
	gdm::Tensor2x2 operator + (gdm::Tensor2x2 a)
	{
		Tensor2x2 p;
		p.xx()=xx_+a.xx();
		p.yy()=yy_+a.yy();
		p.xy()=xy_+a.xy();
		p.yx()=yx_+a.yx();
		return p;
	}

    
  };
} // namespace gdm

#endif // _tensor_hpp



