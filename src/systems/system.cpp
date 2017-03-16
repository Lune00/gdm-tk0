#include "system.hpp"
#include "talk.hpp"
#include "generalCD.hpp"
#include "biaxial.hpp"
#include "biaxial_dof.hpp"
#include "shear.hpp"
#include "shearP_CD.hpp"
#include "check_friction.hpp"
#include "brazilian.hpp"

System* System::factory(string type)
{
  if      (type == "generalCD") return new generalCD; // todo: change this name
  else if (type == "biaxial")   return new Biaxial;
  else if (type == "biaxial_dof")   return new Biaxial_dof;
  else if (type == "shear")     return new Shear;
  else if (type == "shearP_CD") return new shearP_CD;
  else if (type == "check_friction") return new Check_friction;
  else if (type == "brazilian") return new brazilian;

  cerr << "@System::factory, Unknown type of system: " << type << endl << flush;
  gdm::fatal("No System is defined!");
  return 0;
}
