#include "latticealgo.hpp"

using namespace std;

void LatticeAlgo::read_parameters(istream & is)
{
  string token;
  
  is >> token;	
  while(is)
  {	
    if      (token == "dt")      is >> _dt;
    else if (token == "ns")      is >> _ns;
    else if (token == "nsi")     is >> _nsi;
    else if (token == "nsf")     is >> _nsf;
    else if (token == "nitermn") is >> _nitermn;
    else if (token == "nitermx") is >> _nitermx;
    else if (token == "}")       break;
    else cerr << "Unknown parameter: " << token << endl;
    
    is >> token;
  }
}



// Conjugate Gradient iteration
unsigned int LatticeAlgo::iter()
{
//...
  
  return _nitermx;
}


// Initialize the computation
void LatticeAlgo::stand()
{
	
}
