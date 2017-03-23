#ifndef _speedProfile_hpp
#define _speedProfile_hpp

#include "sample.hpp"
#include "network.hpp"
#include "heightProbe.hpp"
#include <vector>


using namespace std;

void speedProfile( vector <heightProbe*> &,vector <double> &,vector <double> &, Sample&  );

void zProfile( vector <heightProbe*> &,vector <double> &, Sample&  );

#endif // _speedProfile_hpp

