#include "inter2d.hpp"
#include "dkdk.hpp"
#include "dkdkP.hpp"
#include "pgpg.hpp"
//#include "rlrl.hpp"
#include "dkrl.hpp"
#include "dkpg.hpp"
#include "pgrl.hpp"

inter2d* inter2d::factory(body2d* b1, body2d* b2)
{
  // dkdk
	if ((typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(disk)))
		return new dkdk(dynamic_cast<disk*>(b1),dynamic_cast<disk*>(b2));

  // pgpg
	if ((typeid(*b1) == typeid(polyg)) && (typeid(*b2) == typeid(polyg)))
		return new pgpg(dynamic_cast<polyg*>(b1),dynamic_cast<polyg*>(b2));
  
  //dkrl
	if ((typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(rline)))
		return new dkrl(dynamic_cast<disk*>(b1),dynamic_cast<rline*>(b2),false);
	if ((typeid(*b1) == typeid(rline)) && (typeid(*b2) == typeid(disk)))
		return new dkrl(dynamic_cast<disk*>(b2),dynamic_cast<rline*>(b1),true);
		
  //pgrl
	if ((typeid(*b1) == typeid(polyg)) && (typeid(*b2) == typeid(rline)))
		return new pgrl(dynamic_cast<polyg*>(b1),dynamic_cast<rline*>(b2),false);
		
	if ((typeid(*b1) == typeid(rline)) && (typeid(*b2) == typeid(polyg)))
		return new pgrl(dynamic_cast<polyg*>(b2),dynamic_cast<rline*>(b1),true);

  //dkpg
	if ((typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(polyg)))
		return new dkpg(dynamic_cast<disk*>(b1),dynamic_cast<polyg*>(b2),false);
	if ((typeid(*b1) == typeid(polyg)) && (typeid(*b2) == typeid(disk)))
		return new dkpg(dynamic_cast<disk*>(b2),dynamic_cast<polyg*>(b1),true);
		
	// rlrl
	  //if ((typeid(*b1) == typeid(rline)) && (typeid(*b2) == typeid(rline)))
		//	return new rlrl(dynamic_cast<rline*>(b1),dynamic_cast<rline*>(b2));
	
	cerr << "@inter2d::factory, Unknown type of interaction" << endl << flush;
	return 0;
}

inter2d* inter2d::factoryP(body2d* b1, body2d* b2, double P)
{
  // dkdk
	if ((typeid(*b1) == typeid(disk)) && (typeid(*b2) == typeid(disk)))
		return new dkdkP(dynamic_cast<disk*>(b1),dynamic_cast<disk*>(b2),P);

	cerr << "@inter2d::factory, Unknown type of interaction" << endl << flush;
	return 0;
}

inter2d* inter2d::factory(string type)
{
  if(type == "dkdk")  return new dkdk;
  if(type == "dkdkP") return new dkdkP;
  if(type == "dkrl")  return new dkrl;
  if(type == "dkpg")  return new dkpg;
  if(type == "pgrl")  return new pgrl;
  if(type == "pgpg")  return new pgpg;
	
  cerr << "@inter2d::factory, Unknown type of interaction: " << type << endl << flush;
  return new dkdk;
}


