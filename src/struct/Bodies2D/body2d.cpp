#include "body2d.hpp"
#include "disk.hpp"
#include "rline.hpp"
#include "polyg.hpp"

body2d* body2d::factory(string type)
{
  if(type == "disk")   return new disk;
  if(type == "polyg")  return new polyg;
  if(type == "rline")  return new rline;
	
  cerr << "@body2d::factory, Unknown type of body: " << type << endl << flush;
  return new disk;
}


