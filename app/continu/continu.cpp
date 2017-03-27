#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include "simulation.hpp"
#include "system.hpp"
#include "dataSet.hpp"
#include<algorithm>



class Grid{

	private:
		unsigned int nx_ , ny_ ;
		double xmin_, xmax_ , ymin_, ymax_ ;
		double dx_ , dy_ ;
		double * array ; 
	public:
		Grid(){};
		~Grid();
		Grid(unsigned int, unsigned int);
		double	getPoint(int,int);

};


Grid::Grid(unsigned int nx , unsigned int ny)
{
	nx_ = nx ;
	ny_ = ny ;
	array = new double [ nx_ * ny_ ];
	xmin_ = 0. ;
	ymin_ = 0. ;
	xmax_ = 0. ;
	ymax_ = 0. ;
	dx_ = 0. ;
	dy_ = 0. ;

	for(int i = 0 ; i != nx_ ; i++)
	{
		for(int j = 0 ; j!= ny_ ; j++)
		{

			array[i*ny_ + j ] = 0. ;
		}
	}
}
Grid::~Grid()
{
	delete array;
}

double Grid::getPoint(int i, int j)
{
	if( array != NULL) return array[i*ny_ + j];
	else return 6. ;
}

int main (int argc,char **argv)
{
	cout<<"Hello (continu) world !"<<endl;
	cout<<"Let's start."<<endl;
	Grid grid(128,128);
	cerr<<grid.getPoint(12,12)<<endl;
	return 0;
}


