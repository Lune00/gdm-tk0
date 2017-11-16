#ifndef _point_hpp
#define _point_hpp

#include <utility>
#include <map>
#include "body2d.hpp"
#include "dof.hpp"
#include "disk.hpp" 
#include "polyg.hpp" // ?
#include "rline.hpp" // ?


class Point{

	private :
		int i_ , j_ ;
		double x_ , y_ ;
		int id_ ;
		std::map<int,double> particules_; // id particule / distance au point
	public:
		Point(){x_ = 0. ; y_ = 0. ; id_ = 0;}
		~Point(){ particules_.clear();};
		double getX(){return x_;};
		double getY(){return y_;};
		int geti() const {return i_ ;}
		int getj() const {return j_ ;}
		int getid() const {return id_;}
		void setid(int id) {id_ = id ; }
		void seti(int i) {i_ = i ;}
		void setj(int j) {j_ = j ;}
		void setX(double x){x_ = x ;};
		void setY(double y){y_ = y ;};
		void add(int,double);
		std::map<int,double> getparticules() const {return particules_;}
		void clearparticules() { particules_.clear();}
		int getsizeparticules() {return particules_.size();}
		bool operator==(const Point&a) const{
			if( id_ != a.id_ ) return false;
			else
				return true;
		}
		bool operator<(const Point&a) const{
			return id_ < a.id_ ;
		}
};
#endif // _point_hpp
