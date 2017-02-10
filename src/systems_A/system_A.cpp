#include "system_A.hpp"
#include "talk.hpp"
#include "shearP_CD_A.hpp"
#include "biaxial_A.hpp"
#include "brazilian_A.hpp"

System_A* System_A::factory(string type)
{
 
  if      (type == "shearP_CD_A") return new shearP_CD_A;
  else if (type == "biaxial_A")   return new Biaxial_A;
  else if (type == "brazilian_A")   return new brazilian_A;
 
  cerr << "@System_A ::factory, Unknown type of system_A: " << type << endl << flush;
  gdm::fatal("No required System_A is defined!");
  return 0;
}


