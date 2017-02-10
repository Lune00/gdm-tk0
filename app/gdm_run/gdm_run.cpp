#include "simulation.hpp"

int main(int argc,char **argv)
{
	Simulation * mySimu = new Simulation();
	mySimu->read_data(argv[1]);
	mySimu->init();
	mySimu->run();

	return 0;
}

