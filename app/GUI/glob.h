#include "simulation.hpp"
#include "CDalgo.hpp"

Simulation  Simu;
CDalgo * cdalgo;
unsigned int nwiter=100;

string WorkingPath;

unsigned int stopComputation = 0;
bool dataReadOK = false;
bool HistorySave = false;
bool saveNetwork = false;
unsigned int iHist = 0;

