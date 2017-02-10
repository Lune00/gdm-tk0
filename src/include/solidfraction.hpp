#ifndef _solidfraction_hpp
#define _solidfraction_hpp

#include "sample.hpp"
#include "network.hpp"
#include "heightProbe.hpp"
#include "circularProbe.hpp"
#include "rectangularProbe.hpp"


using namespace std;

double solidFraction( circularProbe &, Sample & , Network & );

double solidFraction( rectangularProbe &, Sample & , Network & );


double solidFraction( heightProbe &, Sample & , Network & );


#endif // _solidfraction_hpp

