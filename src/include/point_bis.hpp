#ifndef _point_hpp
#define _point_hpp

#include <map>
#include "body2d.hpp"
#include "dof.hpp"
#include "disk.hpp" 
#include "polyg.hpp" // ?
#include "rline.hpp" // ?


class Point{

	private :
		double x_ , y_ ;
		std::map<int,double> particules_; // id particule / distance au point
	public:
		Point(){x_ = 0. ; y_ = 0. ;}
		~Point(){};
		double getX(){return x_;};
		double getY(){return y_;};
		void setX(double x){x_ = x ;};
		void setY(double y){y_ = y ;};
};

#endif // _point_hpp
