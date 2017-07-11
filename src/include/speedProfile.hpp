#ifndef _speedProfile_hpp
#define _speedProfile_hpp

#include "sample.hpp"
#include "network.hpp"
#include "system.hpp"
#include "heightProbe.hpp"
#include <vector>


using namespace std;

void speedProfile( vector <heightProbe*> &,vector <double> &,vector <double> &, Sample&, vector <double>&,  bool shearrateProfile);

void zProfile( vector <heightProbe*> &,vector <double> &, Sample&  );
void TemperatureProfile( vector < heightProbe* > & lprb,vector <double> & Xprofile,vector <double> & Yprofile, vector <double> & XYprofile, Sample& spl, System* sys_);

void RotKeProfile( vector < heightProbe* > & lprb,vector <double> & ROT, Sample& spl);
#endif // _speedProfile_hpp

