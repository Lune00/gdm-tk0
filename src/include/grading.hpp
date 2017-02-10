#ifndef _grading_hpp
#define _grading_hpp

#include "body2d.hpp"
#include "pointSet.hpp"


using namespace std;

//! \brief 
//! \author C. Voivret

//Evaluation de la courbe granulometrique cumulé sur Nc classes
pointSet gradingCurve( vector< body2d * > bodies, unsigned int Nc); 

#endif // _stress_hpp
