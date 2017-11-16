//! \file near.hpp

#ifndef _near_h
#define _near_h

#include <typeinfo>
#include "vertex.hpp"
#include "body2d.hpp"
#include "disk.hpp"
#include "polyg.hpp"
#include "rline.hpp"

//! \brief Factory function that call other functions according to their type.
//!        Return true if the distance between the bodies is lower than dv.
//! \author V. Richefeu 
bool near(const body2d*, const body2d*, const double dv);

//! \author V. Richefeu
bool near(const disk*, const disk*, const double dv);

//! \author V. Richefeu
bool near(const polyg* p1, const polyg* p2, const double dv);

//bool near(const rline*, const rline*, const double dv);

//! \author V. Richefeu
bool near(const disk*, const rline*, const double dv);

//! \author C. Voivret
bool near(const polyg*, const rline*, const double dv);

//! \author V. Richefeu
bool near(const disk* d, const polyg* p, const double dv);

#endif // _near_h

