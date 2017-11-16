#ifndef _tensor_hpp
#define _tensor_hpp
#include <math.h>
#include <iostream>
#include <algorithm>


using namespace std;

class Vecteur
{
	private :
		double x_ ;
		double y_ ;
	public:
		Vecteur() : x_(0.), y_(0.){}
		Vecteur(double x,double y) : x_(x), y_(y) {}
		double getNorme() { return sqrt ( (x_ * x_ + y_ * y_) ) ; } 
		double getX() {return x_;}
		double getY() {return y_;}
		void setX(double x) {x_ = x ;}
		void setY(double y) {y_ = y ;}
};








#endif // _tensor_hpp
