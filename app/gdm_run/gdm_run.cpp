#include "simulation.hpp"
int main(int argc,char **argv)
{
	Simulation * mySimu = new Simulation();
	mySimu->read_data(argv[1]);
	mySimu->init();
	mySimu->run();
	delete mySimu ;

	return 0;
}

