#include "point_bis.hpp"


using namespace std;

//Ajoute une particule au point : particule p , distance au point d
void Point::add(int id, double d){
	particules_.insert( std::make_pair<int,double> (id,d) );
}
