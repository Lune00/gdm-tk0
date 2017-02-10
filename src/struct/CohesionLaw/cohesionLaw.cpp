#include "cohesionLaw.hpp"
#include "weightedAdhesion.cpp"
#include "simpleAdhesion.cpp"
#include "nocohesion.cpp"
//#include "Fragmentation.hpp"

cohesionLaw * cohesionLaw::factory(string type)
{
  if      (type == "weightedAdhesion")   return new weightedAdhesion;
  else if (type == "simpleAdhesion") return new simpleAdhesion;
 // else if (type == "Fragmentation") return new Fragmentation;
  else if (type == "nocohesion")   return new nocohesion;
	
  cerr << "@cohesionLaw::factory, Unknown type of law: " << type << endl << flush;
return new nocohesion;
}


